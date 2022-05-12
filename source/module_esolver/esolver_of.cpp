#include "./esolver_of.h"

//-----------temporary-------------------------
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_symmetry/symmetry.h"
// #include "../src_pw/vdwd2.h"
// #include "../src_pw/vdwd3.h"
// #include "../src_pw/vdwd2_parameters.h"
// #include "../src_pw/vdwd3_parameters.h"
#include "../src_pw/pw_complement.h"
#include "../src_pw/pw_basis.h"
#include "../src_pw/symmetry_rho.h"
#include "../src_io/print_info.h"
#include "../src_pw/H_Ewald_pw.h"
#include "../src_pw/electrons.h"
// #include "../src_pw/occupy.h"
// #include "../src_io/chi0_standard.h"
// #include "../src_io/chi0_hilbert.h"
// #include "../src_io/epsilon0_pwscf.h"
// #include "../src_io/epsilon0_vasp.h"
//-----force-------------------
#include "../src_pw/forces.h"
//-----stress------------------
#include "../src_pw/stress_pw.h"
//---------------------------------------------------

namespace ModuleESolver
{

void ESolver_OF::Init(Input &inp, UnitCell_pseudo &cell)
{
    // save necessary parameters
    this->of_kinetic = inp.of_kinetic;
    this->of_method = inp.of_method;
    this->of_conv = inp.of_conv;
    this->of_tole = inp.of_tole;
    this->of_tolp = inp.of_tolp;
    // this->nelec = inp.nelec;
    this->maxIter = inp.scf_nmax;

    GlobalC::CHR.cal_nelec();
    this->nelec = GlobalC::CHR.nelec;

	if(GlobalC::ucell.atoms[0].xc_func=="HSE"||GlobalC::ucell.atoms[0].xc_func=="PBE0")
	{
		XC_Functional::set_xc_type("pbe");
	}
	else
	{
		XC_Functional::set_xc_type(GlobalC::ucell.atoms[0].xc_func);
	}

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SETUP UNITCELL");
    
    // symmetry analysis should be performed every time the cell is changed
    if (ModuleSymmetry::Symmetry::symm_flag)
    {
        GlobalC::symm.analy_sys(GlobalC::ucell, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "SYMMETRY");
    }

    // Setup the k points according to symmetry.
    GlobalC::kv.set( GlobalC::symm, GlobalV::global_kpoint_card, GlobalV::NSPIN, GlobalC::ucell.G, GlobalC::ucell.latvec );
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT K-POINTS");

    // print information
    // mohan add 2021-01-30
    Print_Info::setup_parameters(GlobalC::ucell, GlobalC::kv);

    // Initalize the plane wave basis set
    GlobalC::pw.gen_pw(GlobalV::ofs_running, GlobalC::ucell, GlobalC::kv);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running,"INIT PLANEWAVE");
    std::cout << " UNIFORM GRID DIM     : " << GlobalC::pw.nx <<" * " << GlobalC::pw.ny <<" * "<< GlobalC::pw.nz << std::endl;
    std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::pw.nbx <<" * " << GlobalC::pw.nby <<" * "<< GlobalC::pw.nbz << std::endl;

    // mohan add 2010-09-13
    // initialize the real-space uniform grid for FFT and parallel
    // distribution of plane waves
    GlobalC::Pgrid.init(GlobalC::pw.ncx, GlobalC::pw.ncy, GlobalC::pw.ncz, GlobalC::pw.nczp,
    GlobalC::pw.nrxx, GlobalC::pw.nbz, GlobalC::pw.bz); // mohan add 2010-07-22, update 2011-05-04
 
    // basis_pw.init(inp, cell);
    this->nrxx = GlobalC::pw.nrxx;
    this->dV = cell.omega / GlobalC::pw.ncxyz; // volume of one point !!! MAYBE WRONG !!!
    this->tf.set_para(this->nrxx, this->dV);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    delete[] this->pdLdphi;
    delete[] this->pdirect;
    this->pdLdphi = new double[this->nrxx];
    this->pdirect = new double[this->nrxx];
    ModuleBase::GlobalFunc::ZEROS(this->pdLdphi, this->nrxx);
    ModuleBase::GlobalFunc::ZEROS(this->pdirect, this->nrxx);

    // Calculate Structure factor
    GlobalC::pw.setup_structure_factor();
    // cout<<"after pgrid init nrxx = "<<GlobalC::pw.nrxx<<endl;

    //----------------------------------------------------------
    // 1 read in initial data:
    //   a lattice structure:atom_species,atom_positions,lattice vector
    //   b k_points
    //   c pseudopotential
    // 2 setup planeware basis, FFT,structure factor, ...
    // 3 initialize local and nonlocal pseudopotential in G_space
    // 4 initialize charge desity and warefunctios in G_space
    //----------------------------------------------------------

    //=====================================
    // init charge/potential/wave functions
    //=====================================
    GlobalC::CHR.allocate(GlobalV::NSPIN, GlobalC::pw.nrxx, GlobalC::pw.ngmc);
    GlobalC::pot.allocate(GlobalC::pw.nrxx);

    GlobalC::wf.allocate(GlobalC::kv.nks);
    this->ppsi = psi::Psi<double>(1, GlobalV::NSPIN, this->nrxx);
    // this->ppsi = new psi::Psi<double>(1, GlobalV::NSPIN, this->nrxx);
    // this->ppsi = new psi::Psi<double>(GlobalC::pw);
    // this->ppsi->resize(1, GlobalV::NSPIN, this->nrxx);
 
    // cout<<GlobalC::pw.nrxx<<endl;
    // cout<<"before ufft allocate"<<endl;
    GlobalC::UFFT.allocate();

    // cout<<"after ufft allocate"<<endl;

    //=======================
    // init pseudopotential
    //=======================
    GlobalC::ppcell.init(GlobalC::ucell.ntype);

    //=====================
    // init hamiltonian
    // only allocate in the beginning of ELEC LOOP!
    //=====================
    GlobalC::hm.hpw.allocate(GlobalC::wf.npwx, GlobalV::NPOL, GlobalC::ppcell.nkb, GlobalC::pw.nrxx);

    //=================================
    // initalize local pseudopotential
    //=================================
    GlobalC::ppcell.init_vloc(GlobalC::pw.nggm, GlobalC::ppcell.vloc);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LOCAL POTENTIAL");

    //======================================
    // Initalize non local pseudopotential
    //======================================
    GlobalC::ppcell.init_vnl(GlobalC::ucell);
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "NON-LOCAL POTENTIAL");

    //=========================================================
    // calculate the total local pseudopotential in real space
    //=========================================================
    GlobalC::pot.init_pot(0, GlobalC::pw.strucFac); //atomic_rho, v_of_rho, set_vrs

    GlobalC::pot.newd();

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT POTENTIAL");

    //==================================================
    // create GlobalC::ppcell.tab_at , for trial wave functions.
    //==================================================
    // GlobalC::wf.init_at_1();

    //================================
    // Initial start wave functions
    //================================
    // if (GlobalV::NBANDS != 0 || (GlobalV::CALCULATION != "scf-sto" && GlobalV::CALCULATION != "relax-sto" && GlobalV::CALCULATION != "md-sto")) //qianrui add
    // {
    //     GlobalC::wf.wfcinit();
    // }
    //===================================No SPIN yet======================================
    for (int ibs = 0; ibs < this->nrxx; ++ibs)
    {
        this->ppsi(0, 0, ibs) = sqrt(GlobalC::CHR.rho[0][ibs]);
    }
    delete[] this->pphi;
    this->pphi = this->ppsi.get_pointer(0);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    // for (int i = 0 ; i < GlobalC::pw.nrxx; ++i)
    // {
    //     cout << GlobalC::CHR.rho[0][i];
    // }
    // p_es->init(inp, cell, basis_pw);
    // p_sqrtRho = sqrt(pes.rho); // pseudo code

    // phamilt->init(bas); 
    // phamilt->initpot(cell, pes);
    if (this->of_method == "tn")
    {
        this->opt_tn.allocate(this->nrxx);
    }
    else if (this->of_method == "cg1" || this->of_method == "cg2")
    {
        this->opt_cg.allocate(this->nrxx);
    }
    else if (this->of_method == "bfgs")
    {
        return;
    }
}

void ESolver_OF::Run(int istep, UnitCell_pseudo& cell)
{
    // =================== for test ===================
    cout << "ncxyz" << setw(4) << GlobalC::pw.ncxyz << endl;
    cout << "nrxx" << setw(4) << this->nrxx << endl;
    cout << "dV" << setw(4) << this->dV << endl;
    cout << "kedf" << setw(4) << this->of_kinetic << endl;
    cout << "nelec" << setw(4) << this->nelec << endl;

    H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::pw);
    // Symmetry_rho srho;
    // for(int is=0; is<GlobalV::NSPIN; is++)
    // {
    //     srho.begin(is, GlobalC::CHR,GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
    // }
    cout << "============== OFDFT BEGIN ==============" <<  endl;
    while(true)
    {
        this->updateV();
        cout << "============== UPDATEV DONE ==============" <<  endl;
        cout << "chemical potential    " << this->mu << endl;
        this->cal_Energy(GlobalC::en);
        this->energy_llast = this->energy_last;
        this->energy_last = this->energy_current;
        this->energy_current = GlobalC::en.etot;
        cout << "TOTAL ENERGY" << "\t" << energy_current << endl;
        if (this->checkExit())
        {
            break;
        }
        cout << "============== CHECK DONE ==============" <<  endl;
        this->solveV();
        cout << "============== SOLVEV DONE ==============" <<  endl;
        this->updateRho();
        cout << "============== UPDATE DONE ==============" <<  endl;
        this->iter++;
        break;
    }
}

// 
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * phi
// 
void ESolver_OF::updateV()
{
    double *dEdphi = new double[this->nrxx];
    // =============================================
    // Be careful !!
    // =============================================
    GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    tf.get_potential(GlobalC::CHR.rho);
    for (int i = 0; i < this->nrxx; ++i)
    {
        
        dEdphi[i] = (GlobalC::pot.vr_eff(0,i) + tf.potential[i]) * 2. * this->ppsi(0,0,i);
    }

    this->mu = this->cal_mu(this->pphi, dEdphi);
    for (int i = 0; i < this->nrxx; ++i)
    {
        this->pdLdphi[i] = dEdphi[i] - 2. * this->mu * this->ppsi(0,0,i);
    }

    // ============= for test ==================
    cout << "VEFF   VTF     RHO     DLDPHI  DEDPHI  PHI PSI" << endl;
    for (int i = 0; i < this->nrxx; i+=100)
    {
        cout << GlobalC::pot.vr_eff(0, i) << "\t"
        << tf.potential[i] << setw(12) << "\t"
        << GlobalC::CHR.rho[0][i] << "\t"
        << pdLdphi[i] << "\t"
        << dEdphi[i] << "\t"
        << this->pphi[i] << "\t"
        << this->ppsi(0,0,i) << endl;
    }
    delete[] dEdphi;
}

//
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * ptempPhi and store it in rdLdphi
// 
void ESolver_OF::calV(double *ptempPhi, double *rdLdphi)
{
    double *dEdtempPhi = new double[this->nrxx];
    double **tempRho = new double*[1];
    tempRho[0] = new double[this->nrxx];
    for (int i = 0; i < this->nrxx; ++i)
    {
        tempRho[0][i] = ptempPhi[i] * ptempPhi[i];
    }
    // =============================================
    // Be careful !!
    // =============================================
    GlobalC::pot.vr = GlobalC::pot.v_of_rho(tempRho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    tf.get_potential(tempRho);
    for (int i = 0; i < this->nrxx; ++i)
    {
        dEdtempPhi[i] = (GlobalC::pot.vr_eff(0,i) + tf.potential[i]) * 2. * ptempPhi[i];
    }
    double tempMu = this->cal_mu(ptempPhi, dEdtempPhi);
    for (int i = 0; i < this->nrxx; ++i)
    {
        rdLdphi[i] = dEdtempPhi[i] - 2. * tempMu * ptempPhi[i];
    }
    delete[] dEdtempPhi;
    delete[] tempRho[0];
    delete[] tempRho;
} 

//
// Get optimization direction d and step theta
//
void ESolver_OF::solveV()
{
    // delete[] this->pphi;
    // this->pphi = this->ppsi.get_pointer(0); // spin degenerate

    // (1) get |d0> with optimization algorithm
    if (this->of_method == "tn")
    {
        opt_tn.next_direct(this->pphi, this->pdLdphi, flag, this->pdirect, this, &ESolver_OF::calV);
    }
    else if (this->of_method == "cg1")
    {
        opt_cg.next_direct(this->pphi, 1, this->pdirect);
    }
    else if (this->of_method == "cg2")
    {
        opt_cg.next_direct(this->pphi, 2, this->pdirect);
    }
    else if (this->of_method == "bfgs")
    {
        return;
    }
    else
    {
        ModuleBase::WARNING_QUIT("ESolver_OF", "method must be CG, TN, or BFGS.");
    }
    cout << "==================== TN DONE =========================" << endl;
    for (int i = 0; i < this->nrxx; i+=100)
    {
        cout << this->pdirect[i] << endl;
    }

    double tempTheta = 0; // tempTheta = |d'|/|d0 + phi|, theta = min(theta, tempTheta)

    // (2) make direction orthogonal to phi
    // |d'> = |d0> - |phi><phi|d0>/nelec
    double innerPhiDir = this->inner_product(this->pdirect, this->pphi, this->nrxx, dV=this->dV);
    for (int i = 0; i < this->nrxx; ++i)
    {
        tempTheta += pow(this->pdirect[i] + this->pphi[i], 2);
        this->pdirect[i] = this->pdirect[i] - this->pphi[i] * innerPhiDir / this->nelec;
    }
    tempTheta = 1. / tempTheta;

    // (3) renormalize direction
    // |d> = |d'> * \sqrt(nelec) / <d'|d'>
    double normDir = sqrt(this->inner_product(this->pdirect, this->pdirect, this->nrxx, dV=this->dV));
    for (int i = 0; i < this->nrxx; ++i)
    {
        tempTheta += tempTheta * this->pdirect[i] * this->pdirect[i];
        this->pdirect[i] = sqrt(this->nelec) * this->pdirect[i] / normDir;
    }
    tempTheta = sqrt(tempTheta);

    // (4) line search to find best theta
    this->theta = min(this->theta, tempTheta);
    double *ptempPhi = new double[this->nrxx];
    double **ptempRho = new double*[1];
    ptempRho[0] = new double[this->nrxx]; // No SPIN
    for (int i = 0; i < this->nrxx; ++i)
    {
        ptempPhi[i] = this->pphi[i] * cos(this->theta) + this->pdirect[i] * sin(this->theta);
        ptempRho[0][i] = ptempPhi[i] * ptempPhi[i];
    }
    double E = 0.;
    double dEdtheta = 0.;
    this->task[0] = 'S'; this->task[1] = 'T'; this->task[2] = 'A'; this->task[3] = 'R'; this->task[4] = 'T';
    while (true)
    {
        dEdtheta = this->caldEdtheta(ptempPhi, ptempRho);
        GlobalC::en.calculate_etot();
        E = GlobalC::en.etot;
        E += this->tf.get_energy(ptempRho);
        this->opt_dcsrch.dcSrch(&E, &dEdtheta, &(this->theta), this->task);
        if (this->task[0] == 'F' && this->task[1] == 'G')
        {
            for (int i = 0; i < this->nrxx; ++i)
            {
                ptempPhi[i] = this->pphi[i] * cos(this->theta) + this->pdirect[i] * sin(this->theta);
                ptempRho[0][i] = ptempPhi[i] * ptempPhi[i];
            }
        }
        else if (task[0] == 'C' && task[1] == 'O')
        {
            break;
        }
        else if (task[0] == 'W' && task[1] == 'A')
        {
            GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << task << std::endl; 
            cout << task;
            break;
        } 
        else if (task[0] == 'E' && task[1] == 'R')
        {
            GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << task << std::endl; 
            break;
        }
    }
    delete[] ptempPhi;
    delete[] ptempRho[0];
    delete[] ptempRho;
}

void ESolver_OF::updateRho()
{
    for (int i = 0; i < this->nrxx; ++i)
    {
        this->ppsi(0,0,i) = this->ppsi(0,0,i) * cos(this->theta) + this->pdirect[i] * sin(this->theta);
        GlobalC::CHR.rho[0][i] = this->ppsi(0,0,i) * this->ppsi(0,0,i); // not sure
    }
}

// 
// Calculate chemical potential mu.
// mu = <dE/dphi|phi> / 2nelec.
// 
double ESolver_OF::cal_mu(double *pphi, double *pdEdphi)
{
    double mu = this->inner_product(pphi, pdEdphi, this->nrxx, this->dV);
    cout << "ORGIN MU" << mu << endl;
    mu = mu / (2.0*this->nelec);
    return mu;
}

void ESolver_OF::cal_Energy(energy &en)
{
    en.calculate_etot();
    cout << "ELOC" << "\t" << en.etot << endl;
    double eTF = tf.get_energy(GlobalC::CHR.rho);
    en.etot += eTF;
    cout << "ETF" << "\t" << eTF << endl;
}

//
// Check convergence, return ture if converge or iter >= maxIter
// 
bool ESolver_OF::checkExit()
{
    if (this->iter >= this->maxIter)
    {
        return true;
    }
    
    bool potConv = false;
    bool energyConv = false;

    double normdLdphi = sqrt(this->inner_product(this->pdLdphi, this->pdLdphi, this->nrxx));
    if (normdLdphi < this->of_tolp)
        potConv = true;

    if (this->iter >= 3
        && abs(this->energy_current - this->energy_last) < this->of_tole
        && abs(this->energy_current - this->energy_llast) < this->of_tole)
        energyConv = true;

    if (this->of_conv == "enery" && energyConv)
    {
        this->conv = true;
        return true;
    }
    else if (this->of_conv == "potential" && potConv)
    {
        this->conv = true;
        return true;
    }
    else if (this->of_conv == "both" && potConv && energyConv)
    {
        this->conv = true;
        return true;
    }
    else
    {
        return false;
    }
}

//
// Calculate dE/dTheta
// dE/dTheta = <dE/dtempPhi|dtempPhi/dTheta>
//           = <dE/dtempPhi|-phi*sin(theta)+d*cos(theta)>
//
double ESolver_OF::caldEdtheta(double *ptempPhi, double **ptempRho)
{
    double *pdEdtempPhi = new double[this->nrxx];
    double *pdPhidTheta = new double[this->nrxx];

    // =============================================
    // Be careful !!
    // =============================================
    GlobalC::pot.vr = GlobalC::pot.v_of_rho(ptempRho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    tf.get_potential(ptempRho);
    for (int i = 0; i < this->nrxx; ++i)
    {
        pdEdtempPhi[i] = (GlobalC::pot.vr_eff(0,i) + tf.potential[i]) * 2 * ptempPhi[i];
    }
    for (int i = 0; i < this->nrxx; ++i)
    {
        pdPhidTheta[i] = - this->pphi[i] * sin(this->theta) + this->pdirect[i] * cos(this->theta);
    }
    double dEdtheta = this->inner_product(pdEdtempPhi, pdPhidTheta, this->nrxx, this->dV);
    delete[] pdEdtempPhi;
    delete[] pdPhidTheta;

    return dEdtheta;
}
}