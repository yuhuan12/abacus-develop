#include "module_hsolver/diago_cg.h"

#include "diago_iter_assist.h"
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_base/memory.h"

using namespace hsolver;

template<typename FPTYPE, typename Device>
DiagoCG<FPTYPE, Device>::DiagoCG(const FPTYPE* precondition_in)
{
    this->device = psi::device::get_device_type<Device>(this->ctx);
    this->precondition = precondition_in;
    test_cg = 0;
    reorder = false;
    this->one = new std::complex<FPTYPE>(1.0, 0.0);
    this->zero = new std::complex<FPTYPE>(0.0, 0.0);
    this->neg_one = new std::complex<FPTYPE>(-1.0, 0.0);
}

template<typename FPTYPE, typename Device>
DiagoCG<FPTYPE, Device>::~DiagoCG() {
    // delete this->cg;
    // delete this->phi_m;
    delmem_complex_op()(this->ctx, this->sphi);
    delmem_complex_op()(this->ctx, this->hphi);
    delmem_complex_op()(this->ctx, this->scg);
    delmem_complex_op()(this->ctx, this->pphi);
    delmem_complex_op()(this->ctx, this->gradient);
    delmem_complex_op()(this->ctx, this->g0);
    delmem_complex_op()(this->ctx, this->lagrange);
    delete this->one;
    delete this->zero;
    delete this->neg_one;
}

template<typename FPTYPE, typename Device>
void DiagoCG<FPTYPE, Device>::diag_mock(hamilt::Hamilt<FPTYPE, Device> *phm_in, psi::Psi<std::complex<FPTYPE>, Device> &phi, FPTYPE *eigenvalue_in)
{
    ModuleBase::TITLE("DiagoCG", "diag_once");
    ModuleBase::timer::tick("DiagoCG", "diag_once");
    
    /// out : record for states of convergence
    this->notconv = 0;

    /// initialize variables
    this->dim = phi.get_current_nbas();
    this->dmx = phi.get_nbasis();
    this->n_band = phi.get_nbands();
    this->eigenvalue = eigenvalue_in;

    setmem_var_h_op()(this->cpu_ctx, this->eigenvalue, 0, this->n_band);
    // ModuleBase::GlobalFunc::ZEROS(this->eigenvalue, this->n_band);

    /// record for how many loops in cg convergence
    FPTYPE avg = 0.0;

    //-------------------------------------------------------------------
    // "poor man" iterative diagonalization of a complex hermitian matrix
    // through preconditioned conjugate gradient algorithm
    // Band-by-band algorithm with minimal use of memory
    // Calls hPhi and sPhi to calculate H|phi> and S|phi>
    // Works for generalized eigenvalue problem (US pseudopotentials) as well
    //-------------------------------------------------------------------
    this->phi_m = new psi::Psi<std::complex<FPTYPE>, Device>(phi, 1, 1);
    // this->hphi.resize(this->dmx, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->hphi, this->dmx);
    setmem_complex_op()(this->ctx, this->hphi, 0, this->dmx);
    // this->sphi.resize(this->dmx, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->sphi, this->dmx);
    setmem_complex_op()(this->ctx, this->sphi, 0, this->dmx);

    this->cg = new psi::Psi<std::complex<FPTYPE>, Device>(phi, 1, 1);
    // this->scg.resize(this->dmx, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->scg, this->dmx);
    setmem_complex_op()(this->ctx, this->scg, 0, this->dmx);
    // this->pphi.resize(this->dmx, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->pphi, this->dmx);
    setmem_complex_op()(this->ctx, this->pphi, 0, this->dmx);

    //in band_by_band CG method, only the first band in phi_m would be calculated
    psi::Range cg_hpsi_range(0);

    // this->gradient.resize(this->dmx, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->gradient, this->dmx);
    setmem_complex_op()(this->ctx, this->gradient, 0, this->dmx);
    // this->g0.resize(this->dmx, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->g0, this->dmx);
    setmem_complex_op()(this->ctx, this->g0, 0, this->dmx);
    // this->lagrange.resize(this->n_band, ModuleBase::ZERO);
    resmem_complex_op()(this->ctx, this->lagrange, this->n_band);
    setmem_complex_op()(this->ctx, this->lagrange, 0, this->n_band);

    ModuleBase::Memory::record("DiagoCG", this->dmx * 8);

    for (int m = 0; m < this->n_band; m++)
    {
        if (test_cg > 2)
            GlobalV::ofs_running << "Diagonal Band : " << m << std::endl;
        //copy psi_in into internal psi, m=0 has been done in Constructor
        if(m>0)
        {
            const std::complex<FPTYPE>* psi_m_in = &(phi(m, 0));
            auto pphi_m = this->phi_m->get_pointer();
            // haozhihan replace COPY_ARRAY
            syncmem_complex_op()(this->ctx, this->ctx, pphi_m, psi_m_in, this->dim);
            // ModuleBase::GlobalFunc::COPYARRAY(psi_m_in, pphi_m, this->dim);
        }
        phm_in->sPsi(this->phi_m->get_pointer(), this->sphi, static_cast<size_t>(this->dim)); // sphi = S|psi(m)>
        this->schmit_orth(m, phi);
        phm_in->sPsi(this->phi_m->get_pointer(), this->sphi, static_cast<size_t>(this->dim)); // sphi = S|psi(m)>

        //do hPsi, actually the result of hpsi stored in Operator,
        //the necessary of copying operation should be checked later
        hpsi_info cg_hpsi_in(this->phi_m, cg_hpsi_range, this->hphi);
        phm_in->ops->hPsi(cg_hpsi_in);

        this->eigenvalue[m] = 
            zdot_real_op()(
                this->ctx, 
                this->dim, 
                this->phi_m->get_pointer(), 
                this->hphi);

        int iter = 0;
        FPTYPE gg_last = 0.0;
        FPTYPE cg_norm = 0.0;
        FPTYPE theta = 0.0;
        bool converged = false;
        for (iter = 0; iter < DiagoIterAssist<FPTYPE, Device>::PW_DIAG_NMAX; iter++)
        {
            this->calculate_gradient();
            this->orthogonal_gradient(phm_in, phi, m);
            this->calculate_gamma_cg(iter, gg_last, cg_norm, theta);
            
            hpsi_info cg_hpsi_in(this->cg, cg_hpsi_range, this->pphi);
            phm_in->ops->hPsi(cg_hpsi_in);

            phm_in->sPsi(this->cg->get_pointer(), this->scg, (size_t)this->dim);
            converged = this->update_psi(cg_norm, theta, this->eigenvalue[m]);

            if (converged) {
                break;
            }
        } // end iter

        std::complex<FPTYPE>* psi_temp = &(phi(m, 0));
        // ModuleBase::GlobalFunc::COPYARRAY(this->phi_m->get_pointer(), psi_temp, this->dim);
        syncmem_complex_op()(this->ctx, this->ctx, psi_temp, this->phi_m->get_pointer(), this->dim);

        if (!converged)
        {
            ++this->notconv;
        }
        avg += static_cast<FPTYPE>(iter) + 1.00;

        // reorder eigenvalues if they are not in the right order
        // (this CAN and WILL happen in not-so-special cases)

        if (m > 0 && reorder)
        {
            ModuleBase::GlobalFunc::NOTE("reorder bands!");
            if (eigenvalue[m] - eigenvalue[m - 1] < -2.0 * DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR)
            {
                // if the last calculated eigenvalue is not the largest...
                int i = 0;
                for (i = m - 2; i >= 0; i--)
                {
                    if (eigenvalue[m] - eigenvalue[i] > 2.0 * DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR)
                        break;
                }
                i++;

                // last calculated eigenvalue should be in the i-th position: reorder
                FPTYPE e0 = eigenvalue[m];
                // ModuleBase::GlobalFunc::COPYARRAY(psi_temp, pphi, this->dim);
                syncmem_complex_op()(this->ctx, this->ctx, pphi, psi_temp, this->dim);

                for (int j = m; j >= i + 1; j--)
                {
                    eigenvalue[j] = eigenvalue[j - 1];
                    std::complex<FPTYPE>* phi_j = &phi(j, 0);
                    std::complex<FPTYPE>* phi_j1 = &phi(j-1, 0);
                    // ModuleBase::GlobalFunc::COPYARRAY(phi_j1, phi_j, this->dim);
                    syncmem_complex_op()(this->ctx, this->ctx, phi_j, phi_j1, this->dim);
                }

                eigenvalue[i] = e0;
                std::complex<FPTYPE>* phi_pointer = &phi(i, 0);
                // ModuleBase::GlobalFunc::COPYARRAY(pphi, phi_pointer, this->dim);
                syncmem_complex_op()(this->ctx, this->ctx, phi_pointer, pphi, this->dim);
            } // endif
        } // end reorder
    } // end m

    avg /= this->n_band;
    DiagoIterAssist<FPTYPE, Device>::avg_iter += avg;

    delete this->phi_m;
    delete this->cg;

    ModuleBase::timer::tick("DiagoCG", "diag_once");
} // end subroutine ccgdiagg

template<typename FPTYPE, typename Device>
void DiagoCG<FPTYPE, Device>::calculate_gradient()
{
    if (this->test_cg == 1) {
        ModuleBase::TITLE("DiagoCG", "calculate_gradient");
    }
    // for (int i = 0; i < this->dim; i++)
    // {
    //     //(2) PH|psi>
    //     this->gradient[i] = this->hphi[i] / this->precondition[i];
    //     //(3) PS|psi>
    //     this->pphi[i] = this->sphi[i] / this->precondition[i];
    // }
    // denghui replace this at 20221106
    // TODO: use GPU precondition to initialize CG class
    if (this->device == psi::GpuDevice) {
        vector_div_vector_op<FPTYPE, Device>()(this->ctx, this->dim, this->gradient, this->hphi, this->d_precondition);
        vector_div_vector_op<FPTYPE, Device>()(this->ctx, this->dim, this->pphi, this->sphi, this->d_precondition);
    }
    else {
        vector_div_vector_op<FPTYPE, Device>()(this->ctx, this->dim, this->gradient, this->hphi, this->precondition);
        vector_div_vector_op<FPTYPE, Device>()(this->ctx, this->dim, this->pphi, this->sphi, this->precondition);
    }

    // Update lambda !
    // (4) <psi|SPH|psi >
    const FPTYPE eh = hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, this->sphi, this->gradient);
    // (5) <psi|SPS|psi >
    const FPTYPE es = hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, this->sphi, this->pphi);
    const FPTYPE lambda = eh / es;

    // Update g!
    // for (int i = 0; i < this->dim; i++)
    // {
    //     //               <psi|SPH|psi>
    //     // (6) PH|psi> - ------------- * PS |psi>
    //     //               <psi|SPS|psi>
    //     //
    //     // So here we get the gradient.
    //     this->gradient[i] -= lambda * this->pphi[i];
    // }
    // haozhihan replace this 2022-10-6
    constantvector_addORsub_constantVector_op<FPTYPE, Device>()(this->ctx, this->dim, this->gradient, this->gradient, 1.0, this->pphi, (-lambda));
}

template<typename FPTYPE, typename Device>
void DiagoCG<FPTYPE, Device>::orthogonal_gradient(hamilt::Hamilt<FPTYPE, Device> *phm_in, const psi::Psi<std::complex<FPTYPE>, Device> &eigenfunction, const int m)
{
    if (test_cg == 1) {
        ModuleBase::TITLE("DiagoCG", "orthogonal_gradient");
    }
    // ModuleBase::timer::tick("DiagoCG","orth_grad");

    phm_in->sPsi(this->gradient, this->scg, (size_t)this->dim);
    // int inc = 1;
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-07
    gemv_op<FPTYPE, Device>()(
        this->ctx,
        'C',
        this->dim,
        m,
        this->one,
        eigenfunction.get_pointer(),
        this->dmx,
        this->scg,
        1,
        this->zero,
        this->lagrange,
        1);

    Parallel_Reduce::reduce_complex_double_pool(this->lagrange, m);

    // (3) orthogonal |g> and |scg> to all states (0~m-1)
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-07
    gemv_op<FPTYPE, Device>()(
        this->ctx,
        'N',
        this->dim,
        m,
        this->neg_one,
        eigenfunction.get_pointer(),
        this->dmx,
        this->lagrange,
        1,
        this->one,
        this->gradient,
        1);

    gemv_op<FPTYPE, Device>()(
        this->ctx,
        'N',
        this->dim,
        m,
        this->neg_one,
        eigenfunction.get_pointer(),
        this->dmx,
        this->lagrange,
        1,
        this->one,
        this->scg,
        1);
}

template<typename FPTYPE, typename Device>
void DiagoCG<FPTYPE, Device>::calculate_gamma_cg(const int iter, FPTYPE &gg_last, const FPTYPE &cg_norm, const FPTYPE &theta)
{
    if (test_cg == 1) {
        ModuleBase::TITLE("DiagoCG", "calculate_gamma_cg");
    }
    auto pcg = this->cg->get_pointer();
    auto pphi_m = this->phi_m->get_pointer();
    FPTYPE gg_inter;
    if (iter > 0)
    {
        // (1) Update gg_inter!
        // gg_inter = <g|g0>
        // Attention : the 'g' in g0 is getted last time
        gg_inter
            = hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, this->gradient, this->g0); // b means before
    }

    // (2) Update for g0!
    // two usage:
    // firstly, for now, calculate: gg_now
    // secondly, prepare for the next iteration: gg_inter
    // |g0> = P | scg >
    // for (int i = 0; i < this->dim; i++)
    // {
    //     this->g0[i] = this->precondition[i] * this->scg[i];
    // }
    // denghui replace this 20221106
    // TODO: use GPU precondition instead
    if (this->device == psi::GpuDevice) {
        vector_mul_vector_op<FPTYPE, Device>()(this->ctx, this->dim, this->g0, this->scg, this->d_precondition);
    }
    else {
        vector_mul_vector_op<FPTYPE, Device>()(this->ctx, this->dim, this->g0, this->scg, this->precondition);
    }

    // (3) Update gg_now!
    // gg_now = < g|P|scg > = < g|g0 >
    const FPTYPE gg_now = hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, this->gradient, this->g0);

    if (iter == 0)
    {
        // (40) gg_last first value : equal gg_now
        gg_last = gg_now;
        // (50) cg direction first value : |g>
        // |cg> = |g>
        // ModuleBase::GlobalFunc::COPYARRAY(this->gradient, pcg, this->dim);
        syncmem_complex_op ()(this->ctx, this->ctx, pcg, this->gradient, this->dim);
    }
    else
    {
        // (4) Update gamma !
        assert(gg_last != 0.0);
        const FPTYPE gamma = (gg_now - gg_inter) / gg_last;

        // (5) Update gg_last !
        gg_last = gg_now;

        // (6) Update cg direction !(need gamma and |go> ):
        // for (int i = 0; i < this->dim; i++)
        // {
        //     pcg[i] = gamma * pcg[i] + this->gradient[i];
        // }
        // haozhihan replace this 2022-10-6
        constantvector_addORsub_constantVector_op<FPTYPE, Device>()(this->ctx, this->dim, pcg, pcg, gamma, this->gradient, 1.0);

        const FPTYPE norma = gamma * cg_norm * sin(theta);
        std::complex<FPTYPE> znorma(norma * -1, 0.0);

        // haozhihan replace this 2022-10-6
        // const int one = 1;
        // zaxpy_(&this->dim, &znorma, pphi_m, &one, pcg, &one);
        /*for (int i = 0; i < this->dim; i++)
        {
            pcg[i] -= norma * pphi_m[i];
        }*/
        axpy_op<FPTYPE, Device>()(this->ctx, this->dim, &znorma, pphi_m, 1, pcg, 1);
    }
}

template<typename FPTYPE, typename Device>
bool DiagoCG<FPTYPE, Device>::update_psi(FPTYPE &cg_norm, FPTYPE &theta, FPTYPE &eigenvalue)
{
    if (test_cg == 1)
        ModuleBase::TITLE("DiagoCG", "update_psi");
    cg_norm = sqrt(hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, this->cg->get_pointer(), this->scg));

    if (cg_norm < 1.0e-10)
        return 1;

    std::complex<FPTYPE>* phi_m_pointer = this->phi_m->get_pointer();

    const FPTYPE a0
        = hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, phi_m_pointer, this->pphi) * 2.0 / cg_norm;
    const FPTYPE b0
        = hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, this->dim, this->cg->get_pointer(), this->pphi) / (cg_norm * cg_norm);

    const FPTYPE e0 = eigenvalue;
    theta = atan(a0 / (e0 - b0)) / 2.0;

    const FPTYPE new_e = (e0 - b0) * cos(2.0 * theta) + a0 * sin(2.0 * theta);

    const FPTYPE e1 = (e0 + b0 + new_e) / 2.0;
    const FPTYPE e2 = (e0 + b0 - new_e) / 2.0;

    if (e1 > e2)
    {
        theta += ModuleBase::PI_HALF;
    }

    eigenvalue = min(e1, e2);
    //	OUT("eigenvalue",eigenvalue);

    const FPTYPE cost = cos(theta);
    const FPTYPE sint_norm = sin(theta) / cg_norm;

    auto pcg = this->cg->get_pointer();
    // for (int i = 0; i < this->dim; i++)
    // {
    //     phi_m_pointer[i] = phi_m_pointer[i] * cost + sint_norm * pcg[i];
    // }
    
    // haozhihan replace this 2022-10-6
    constantvector_addORsub_constantVector_op<FPTYPE, Device>()(this->ctx, this->dim, phi_m_pointer, phi_m_pointer, cost, pcg, sint_norm);


    //	std::cout << "\n overlap2 = "  << this->ddot(dim, phi_m, phi_m);

    if (abs(eigenvalue - e0) < DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR)
    {
        // ModuleBase::timer::tick("DiagoCG","update");
        return 1;
    }
    else
    {
        // for (int i = 0; i < this->dim; i++)
        // {
        //     this->sphi[i] = this->sphi[i] * cost + sint_norm * this->scg[i];
        //     this->hphi[i] = this->hphi[i] * cost + sint_norm * this->pphi[i];
        // }

        // haozhihan replace this 2022-10-6
        constantvector_addORsub_constantVector_op<FPTYPE, Device>()(this->ctx, this->dim, this->sphi, this->sphi, cost, this->scg, sint_norm);
        constantvector_addORsub_constantVector_op<FPTYPE, Device>()(this->ctx, this->dim, this->hphi, this->hphi, cost, this->pphi, sint_norm);
        return 0;
    }
}

template<typename FPTYPE, typename Device>
void DiagoCG<FPTYPE, Device>::schmit_orth(
                          const int &m, // end
                          const psi::Psi<std::complex<FPTYPE>, Device> &psi)
{
    //	ModuleBase::TITLE("DiagoCG","schmit_orth");
    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
    // orthogonalize starting eigenfunction to those already calculated
    // phi_m orthogonalize to psi(start) ~ psi(m-1)
    // Attention, the orthogonalize here read as
    // psi(m) -> psi(m) - \sum_{i < m} < psi(i) | S | psi(m) > psi(i)
    // so the orthogonalize is performed about S.
    assert(m >= 0);
    assert(psi.get_nbands() >= m);

    std::complex<FPTYPE> * lagrange_so = nullptr;
    resmem_complex_op()(this->ctx, lagrange_so, m + 1);

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-6
    int inc = 1;
    gemv_op<FPTYPE, Device>()(
        this->ctx,
        'C',
        this->dim,
        (m + 1),
        this->one,
        psi.get_pointer(),
        this->dmx,
        this->sphi,
        inc,
        this->zero,
        lagrange_so,
        inc);

    // be careful , here reduce m+1
    Parallel_Reduce::reduce_complex_double_pool(lagrange_so, m + 1);

    std::complex<FPTYPE> var(0, 0);
    syncmem_complex_d2h_op()(this->cpu_ctx, this->ctx, &var, lagrange_so + m, 1);
    FPTYPE psi_norm = var.real();

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-6
    gemv_op<FPTYPE, Device>()(
        this->ctx,
        'N',
        this->dim,
        m,
        this->neg_one,
        psi.get_pointer(),
        this->dmx,
        lagrange_so,
        inc,
        this->one,
        this->phi_m->get_pointer(),
        inc);

    //======================================================================
    /*for (int j = 0; j < m; j++)
    {
        for (int ig =0; ig < dim; ig++)
        {
            phi_m[ig] -= lagrange[j] * psi(j, ig);
        }
        psi_norm -= ( conj(lagrange[j]) * lagrange[j] ).real();
    }*/
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    psi_norm -= hsolver::zdot_real_op<FPTYPE, Device>()(this->ctx, m, lagrange_so, lagrange_so, false);

    if (psi_norm <= 0.0)
    {
        std::cout << " m = " << m << std::endl;
        for (int j = 0; j <= m; ++j)
        {
            std::cout << "j = " << j << " lagrange norm = " << (conj(lagrange_so[j]) * lagrange_so[j]).real()
                      << std::endl;
        }
        std::cout << " in DiagoCG, psi norm = " << psi_norm << std::endl;
        std::cout << " If you use GNU compiler, it may due to the zdotc is unavailable." << std::endl;
        ModuleBase::WARNING_QUIT("schmit_orth", "psi_norm <= 0.0");
    }

    psi_norm = sqrt(psi_norm);

    auto pphi_m = this->phi_m->get_pointer();
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // haozhihan replace 2022-10-6
    // scal_op<FPTYPE, Device>()(this->ctx, this->dim, &psi_norm, pphi_m, 1);
    //======================================================================
    // for (int ig = 0; ig < this->dim; ig++)
    // {
    //     pphi_m[ig] /= psi_norm;
    // }
    vector_div_constant_op<FPTYPE, Device>()(this->ctx, this->dim, pphi_m, pphi_m, psi_norm);

    // ModuleBase::timer::tick("DiagoCG","schmit_orth");
    delmem_complex_op()(this->ctx, lagrange_so);
}

template<typename FPTYPE, typename Device>
void DiagoCG<FPTYPE, Device>::diag(hamilt::Hamilt<FPTYPE, Device> *phm_in, psi::Psi<std::complex<FPTYPE>, Device> &psi, FPTYPE *eigenvalue_in)
{
    /// record the times of trying iterative diagonalization
    int ntry = 0;
    this->notconv = 0;
    if (this->device == psi::GpuDevice) {
        resmem_var_op()(this->ctx, this->d_precondition, psi.get_nbasis());
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, this->d_precondition, this->precondition, psi.get_nbasis());
    }
    do
    {
        if(DiagoIterAssist<FPTYPE, Device>::need_subspace || ntry > 0)
        {
            DiagoIterAssist<FPTYPE, Device>::diagH_subspace(phm_in, psi, psi, eigenvalue_in);
        }

        DiagoIterAssist<FPTYPE, Device>::avg_iter += 1.0;
        this->reorder = true;

        this->diag_mock(phm_in, psi, eigenvalue_in);

        ++ntry;
    } while (DiagoIterAssist<FPTYPE, Device>::test_exit_cond(ntry, this->notconv));

    if (notconv > max(5, psi.get_nbands() / 4)) {
        std::cout << "\n notconv = " << this->notconv;
        std::cout << "\n DiagoCG::diag', too many bands are not converged! \n";
    }
}

namespace hsolver {
template class DiagoCG<float, psi::DEVICE_CPU>;
template class DiagoCG<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class DiagoCG<float, psi::DEVICE_GPU>;
template class DiagoCG<double, psi::DEVICE_GPU>;
#endif 
} // namespace hsolver