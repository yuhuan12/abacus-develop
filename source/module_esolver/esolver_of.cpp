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
    // this->nelec = GlobalC::CHR.nelec;

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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    // delete[] this->pdLdphi;
    // delete[] this->pdirect;
    this->pdLdphi = new double*[GlobalV::NSPIN];
    this->pdEdphi = new double*[GlobalV::NSPIN];
    this->pdirect = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->pdLdphi[is] = new double[this->nrxx];
        this->pdEdphi[is] = new double[this->nrxx];
        this->pdirect[is] = new double[this->nrxx];
        ModuleBase::GlobalFunc::ZEROS(this->pdLdphi[is], this->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdEdphi[is], this->nrxx);
        ModuleBase::GlobalFunc::ZEROS(this->pdirect[is], this->nrxx);
    }

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
 
    GlobalC::UFFT.allocate();

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
    this->pphi = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ibs = 0; ibs < this->nrxx; ++ibs)
        {
            this->ppsi(0, is, ibs) = sqrt(GlobalC::CHR.rho[is][ibs]);
        }
        this->pphi[is] = this->ppsi.get_pointer(is);
    }

    this->nelec = new double[GlobalV::NSPIN];
    if (GlobalV::NSPIN == 1)
    {
        this->nelec[0] = GlobalC::CHR.nelec;
    }
    else if (GlobalV::NSPIN == 2)
    {
        this->nelec[0] = GlobalC::ucell.magnet.get_nelup();
        this->nelec[1] = GlobalC::ucell.magnet.get_neldw();
    }
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    // p_es->init(inp, cell, basis_pw);
    // p_sqrtRho = sqrt(pes.rho); // pseudo code

    // phamilt->init(bas); 
    // phamilt->initpot(cell, pes);
    if (this->of_method == "tn")
    {
        this->opt_tn.allocate(this->nrxx);
        this->opt_tn.setPara(this->dV);
    }
    else if (this->of_method == "cg1" || this->of_method == "cg2")
    {
        this->opt_cg.allocate(this->nrxx);
        this->opt_cg.setPara(this->dV);
        this->opt_dcsrch.set_paras(1e-4,1e-2);
    }
    else if (this->of_method == "bfgs")
    {
        return;
    }

    if (GlobalV::NSPIN == 2)
    {
        this->opt_cg_mag = new Opt_CG;
        this->opt_cg_mag->allocate(GlobalV::NSPIN);
    }

    this->mu = new double[GlobalV::NSPIN];
    this->theta = new double[GlobalV::NSPIN];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->mu[is] = 0;
        this->theta[is] = 0.;
    }
    if (GlobalV::NSPIN == 1)
    {
        this->theta[0] = 0.2;
    }
}

void ESolver_OF::Run(int istep, UnitCell_pseudo& cell)
{

    H_Ewald_pw::compute_ewald(GlobalC::ucell, GlobalC::pw);
    Symmetry_rho srho;
    for(int is=0; is<GlobalV::NSPIN; is++)
    {
        srho.begin(is, GlobalC::CHR,GlobalC::pw, GlobalC::Pgrid, GlobalC::symm);
    }
    cout << "============================================ OFDFT ============================================" <<  endl;
    while(true)
    {
        this->updateV();
        // cout << "============== UPDATEV DONE ==============" <<  endl;
        this->cal_Energy(GlobalC::en);
        this->energy_llast = this->energy_last;
        this->energy_last = this->energy_current;
        this->energy_current = GlobalC::en.etot;
        this->printInfo();
        if (this->checkExit())
        {
            break;
        }
        this->solveV();
        this->updateRho();
        this->iter++;
    }
    if (GlobalC::CHR.out_chg > 0)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            std::stringstream ssc;
            ssc << GlobalV::global_out_dir << "tmp" << "_SPIN" << is << "_CHG";
            GlobalC::CHR.write_rho(GlobalC::CHR.rho[is], is, iter, ssc.str(), 11);
        }
    }
    if (GlobalV::CAL_FORCE)
    {
        ModuleBase::matrix ff(GlobalC::ucell.nat, 3);
        this->cal_Force(ff);
        // for (int ia = 0; ia < GlobalC::ucell.nat; ++ia)
        // {
        //     for (int i = 0; i < 3; ++i)
        //     {
        //         cout << ff(ia,i) << "    ";
        //     }
        //     cout << endl;
        // }
    }
    if (GlobalV::CAL_STRESS)
    {
        double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8 / 10;
        ModuleBase::matrix stress(3,3);
        this->cal_Stress(stress);
        cout << "STRESS (GPa)" << endl;
        for (int i = 0; i < 3; ++i)
        {
            cout << setprecision(8) << setw(16) << stress(i, 0)*unit_transform << setw(16) << stress(i, 1)*unit_transform << setw(16) << stress(i, 2)*unit_transform << endl;
        }
    }
}

// ===================================================================
// NOTE THIS FUNCTION SHOULD BE UESD AFTER POTENTIAL HAS BEEN UPDATED
// ===================================================================
void ESolver_OF::cal_Energy(energy &en)
{
    en.calculate_etot();
    double eTF = this->tf.get_energy(GlobalC::CHR.rho);
    double ePP = 0.;
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ePP += this->inner_product(GlobalC::pot.vltot, GlobalC::CHR.rho[is], this->nrxx, this->dV);
    }
    en.etot += eTF + ePP;
}

// 
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * phi
// 
void ESolver_OF::updateV()
{
    // double *dEdphi = new double[this->nrxx];

    GlobalC::pot.vr = GlobalC::pot.v_of_rho(GlobalC::CHR.rho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    this->tf.get_potential(GlobalC::CHR.rho);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ModuleBase::GlobalFunc::ZEROS(this->pdEdphi[is], this->nrxx);
        for (int ir = 0; ir < this->nrxx; ++ir)
        { 
            this->pdEdphi[is][ir] = (GlobalC::pot.vr_eff(is,ir) + this->tf.potential[is][ir]) * 2. * this->ppsi(0,is,ir);
        }
        this->mu[is] = this->cal_mu(this->pphi[is], this->pdEdphi[is], this->nelec[is]);

        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pdLdphi[is][ir] = this->pdEdphi[is][ir] - 2. * this->mu[is] * this->ppsi(0,is,ir);
        }
    }

    // delete[] dEdphi;
}

//
// Get optimization direction d and step theta
//
void ESolver_OF::solveV()
{
    // (1) get |d0> with optimization algorithm
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (this->of_method == "tn")
        {
            this->tnSpinFlag = is;
            opt_tn.next_direct(this->pphi[is], this->pdLdphi[is], flag, this->pdirect[is], this, &ESolver_OF::calV);
        }
        else if (this->of_method == "cg1")
        {
            opt_cg.next_direct(this->pdLdphi[is], 1, this->pdirect[is]);
        }
        else if (this->of_method == "cg2")
        {
            opt_cg.next_direct(this->pdLdphi[is], 2, this->pdirect[is]);
        }
        else if (this->of_method == "bfgs")
        {
            return;
        }
        else
        {
            ModuleBase::WARNING_QUIT("ESolver_OF", "method must be CG, TN, or BFGS.");
        }
    }
    // initialize tempPhi and tempRho used in line search
    double **ptempPhi = new double*[GlobalV::NSPIN];
    double **ptempRho = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        ptempPhi[is] = new double[this->nrxx];
        ptempRho[is] = new double[this->nrxx]; // No SPIN
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            ptempPhi[is][ir] = this->pphi[is][ir];
            ptempRho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
        }
    }
    
    // (2) rotate and renormalize the direction
    this->getNextDirect();

    // (3) make sure dEdtheta>0 when theta = 0
    double E = 0.; // energy of tempPhi and tempRho
    double *dEdtheta = new double[GlobalV::NSPIN]; // dE/dtheta at tempPhi
    double *tempTheta = new double[GlobalV::NSPIN];
    ModuleBase::GlobalFunc::ZEROS(dEdtheta, GlobalV::NSPIN);
    ModuleBase::GlobalFunc::ZEROS(tempTheta, GlobalV::NSPIN);
    this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        if (dEdtheta[is] > 0)
        {
            GlobalV::ofs_warning << "ESolver_OF: WARNING " << "dEdphi > 0, replace direct with steepest descent method." << endl;
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                this->pdirect[is][ir] = - this->pdLdphi[is][ir];
            }
            this->getNextDirect();
            this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
            if (dEdtheta[is] > 0)
            {
                GlobalV::ofs_warning << "ESolver_OF: WARNING " << "when use steepes dencent method, dEdphi > 0, so we might get minimum." << endl;
            }
        }
    }
    delete[] tempTheta;

// ======================== for test ============================
    // if (this->iter == 2)
    // {
    //     for (int i = -100; i < 100; ++i)
    //     {
    //         this->theta[0] = 0.001 * i;
    //         for (int i = 0; i < this->nrxx; ++i)
    //         {
    //             ptempPhi[0][i] = this->pphi[0][i] * cos(this->theta[0]) + this->pdirect[0][i] * sin(this->theta[0]);
    //             ptempRho[0][i] = ptempPhi[0][i] * ptempPhi[0][i];
    //         }
    //         this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);
    //         GlobalC::en.calculate_etot();
    //         E = GlobalC::en.etot;
    //         E += this->tf.get_energy(ptempRho);
    //         E += this->inner_product(GlobalC::pot.vltot, ptempRho[0], this->nrxx, this->dV);
    //         GlobalV::ofs_warning << i << "    " << dEdtheta[0] << "    " << E << endl;
    //         if (this->theta[0] == 0) cout << "dEdtheta    " << dEdtheta[0]<< endl;
    //     }
    //     exit(0);
    // }
// ======================== for test ============================

    // (4) line search to find best theta
    if (GlobalV::NSPIN == 1)
    {
        int numDC = 0; // iteration number of line search
        this->task[0] = 'S'; this->task[1] = 'T'; this->task[2] = 'A'; this->task[3] = 'R'; this->task[4] = 'T';
        while (true)
        {
            GlobalC::en.calculate_etot();
            E = GlobalC::en.etot;
            E += this->tf.get_energy(ptempRho);
            E += this->inner_product(GlobalC::pot.vltot, ptempRho[0], this->nrxx, this->dV);
            this->opt_dcsrch.dcSrch(&E, &dEdtheta[0], &(this->theta[0]), this->task);
            numDC++;

            if (this->task[0] == 'F' && this->task[1] == 'G')
            {
                for (int i = 0; i < this->nrxx; ++i)
                {
                    ptempPhi[0][i] = this->pphi[0][i] * cos(this->theta[0]) + this->pdirect[0][i] * sin(this->theta[0]);
                    ptempRho[0][i] = ptempPhi[0][i] * ptempPhi[0][i];
                }
                this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);

                if (numDC > this->maxDCsrch)
                {
                    GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << "excedd the max iter number." << endl;
                    break;
                }
            }
            else if (this->task[0] == 'C' && this->task[1] == 'O')
            {
                break;
            }
            else if (this->task[0] == 'W' && this->task[1] == 'A')
            {
                GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << this->task << std::endl; 
                cout << this->task << endl;
                break;
            } 
            else if (this->task[0] == 'E' && this->task[1] == 'R')
            {
                GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << this->task << std::endl; 
                cout << this->task << endl;
                break;
            }
        }
    }
    else if (GlobalV::NSPIN == 2)
    {
        this->opt_cg_mag->refresh();

        double *pthetaDir = new double[GlobalV::NSPIN];
        double *tempTheta = new double[GlobalV::NSPIN];
        ModuleBase::GlobalFunc::ZEROS(pthetaDir, GlobalV::NSPIN);
        ModuleBase::GlobalFunc::ZEROS(tempTheta, GlobalV::NSPIN);
        double thetaAlpha = 0.;
        double alphaTol = 1e-4;
        double maxThetaDir = 0.;
        double dEdalpha = 0.;
        int thetaIter = 0;
        int numDC = 0;

        while (true)
        {
            this->opt_cg_mag->next_direct(dEdtheta, 1, pthetaDir);

            // this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);
            dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1.);

            if (dEdalpha >= 0.)
            {
                for (int is = 0; is < GlobalV::NSPIN; ++is)
                {
                    pthetaDir[is] = -dEdtheta[is];
                }
                // this->caldEdtheta(ptempPhi, ptempRho, this->theta, dEdtheta);
                dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1);
            }

            maxThetaDir = max(abs(pthetaDir[0]), abs(pthetaDir[1]));
            thetaAlpha = min(0.1, 0.1*ModuleBase::PI/maxThetaDir);

            // line search along thetaDir to find thetaAlpha
            this->opt_dcsrch.set_paras(1e-4, 1e-2, 1e-12, 0., ModuleBase::PI/maxThetaDir);
            this->task[0] = 'S'; this->task[1] = 'T'; this->task[2] = 'A'; this->task[3] = 'R'; this->task[4] = 'T';
            numDC = 0;
            while(true)
            {
                GlobalC::en.calculate_etot();
                E = GlobalC::en.etot;
                E += this->tf.get_energy(ptempRho);
                for (int is = 0; is < GlobalV::NSPIN; ++is) E += this->inner_product(GlobalC::pot.vltot, ptempRho[is], this->nrxx, this->dV);
                this->opt_dcsrch.dcSrch(&E, &dEdalpha, &thetaAlpha, this->task);
                numDC++;

                if (this->task[0] == 'F' && this->task[1] == 'G')
                {
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        tempTheta[is] = this->theta[is] + thetaAlpha * pthetaDir[is];
                        for (int ir = 0; ir < this->nrxx; ++ir)
                        {
                            ptempPhi[is][ir] = this->ppsi(0, is, ir) * cos(tempTheta[is]) + this->pdirect[is][ir] * sin(tempTheta[is]);
                            ptempRho[is][ir] = ptempPhi[is][ir] * ptempPhi[is][ir];
                        }
                    }
                    this->caldEdtheta(ptempPhi, ptempRho, tempTheta, dEdtheta);
                    dEdalpha = this->inner_product(dEdtheta, pthetaDir, 2, 1);

                    if (numDC > 10)
                    {
                        GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << "excedd the max iter number." << endl;
                        break;
                    }
                }
                else if (this->task[0] == 'C' && this->task[1] == 'O')
                {
                    break;
                }
                else if (this->task[0] == 'W' && this->task[1] == 'A')
                {
                    GlobalV::ofs_warning << "ESolver_OF linesearch: WARNING " << this->task << std::endl; 
                    cout << this->task << endl;
                    break;
                } 
                else if (this->task[0] == 'E' && this->task[1] == 'R')
                {
                    GlobalV::ofs_warning << "ESolver_OF linesearch: ERROR " << this->task << std::endl; 
                    cout << this->task << endl;
                    break;
                }
            }

            for (int is = 0; is < GlobalV::NSPIN; ++is) this->theta[is] += thetaAlpha * pthetaDir[is];
            if (sqrt(dEdtheta[0] * dEdtheta[0] + dEdtheta[1] * dEdtheta[1]) < alphaTol) break;
            thetaIter++;
            if (thetaIter > 2) break;
        }
        delete[] tempTheta;
        delete[] pthetaDir;
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] ptempPhi[is];
        delete[] ptempRho[is];
    }
    delete[] ptempPhi;
    delete[] ptempRho;
    delete[] dEdtheta;
}

// 
// rotate and renormalize the direction, make it orthogonal to phi, and <d|d> = nelec
// 
void ESolver_OF::getNextDirect()
{
    if (GlobalV::NSPIN == 1)
    {
        double tempTheta = 0; // tempTheta = |d'|/|d0 + phi|, theta = min(theta, tempTheta)

        // (1) make direction orthogonal to phi
        // |d'> = |d0> - |phi><phi|d0>/nelec
        double innerPhiDir = this->inner_product(this->pdirect[0], this->pphi[0], this->nrxx, dV=this->dV);
        for (int i = 0; i < this->nrxx; ++i)
        {
            tempTheta += pow(this->pdirect[0][i] + this->ppsi(0, 0, i), 2);
            this->pdirect[0][i] = this->pdirect[0][i] - this->ppsi(0, 0, i) * innerPhiDir / this->nelec[0];
        }
        tempTheta = 1. / tempTheta;

        // (2) renormalize direction
        // |d> = |d'> * \sqrt(nelec) / <d'|d'>
        double normDir = sqrt(this->inner_product(this->pdirect[0], this->pdirect[0], this->nrxx, dV=this->dV));
        for (int i = 0; i < this->nrxx; ++i)
        {
            tempTheta += tempTheta * this->pdirect[0][i] * this->pdirect[0][i];
            this->pdirect[0][i] = sqrt(this->nelec[0]) * this->pdirect[0][i] / normDir;
        }
        tempTheta = sqrt(tempTheta);
        this->theta[0] = min(this->theta[0], tempTheta);
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            // (1) make direction orthogonal to phi
            // |d'> = |d0> - |phi><phi|d0>/nelec
            double innerPhiDir = this->inner_product(this->pdirect[is], this->pphi[is], this->nrxx, dV=this->dV);
            for (int i = 0; i < this->nrxx; ++i)
            {
                this->pdirect[is][i] = this->pdirect[is][i] - this->ppsi(0, is, i) * innerPhiDir / this->nelec[is];
            }

            // (2) renormalize direction
            // |d> = |d'> * \sqrt(nelec) / <d'|d'>
            double normDir = sqrt(this->inner_product(this->pdirect[is], this->pdirect[is], this->nrxx, dV=this->dV));
            for (int i = 0; i < this->nrxx; ++i)
            {
                this->pdirect[is][i] = sqrt(this->nelec[is]) * this->pdirect[is][i] / normDir;
            }
            this->theta[is] = 0.;
        }
    }
}

void ESolver_OF::updateRho()
{
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->ppsi(0,is,ir) = this->ppsi(0,is,ir) * cos(this->theta[is]) + this->pdirect[is][ir] * sin(this->theta[is]);
            GlobalC::CHR.rho[is][ir] = this->ppsi(0,is,ir) * this->ppsi(0,is,ir);
        }
    }
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

    this->normdLdphi = 0.;
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
       this->normdLdphi += this->inner_product(this->pdLdphi[is], this->pdLdphi[is], this->nrxx, this->dV);
    }
    this->normdLdphi = sqrt(this->normdLdphi/this->nrxx/GlobalV::NSPIN);

    if (this->normdLdphi < this->of_tolp)
        potConv = true;

    if (this->iter >= 3
        && abs(this->energy_current - this->energy_last) < this->of_tole
        && abs(this->energy_current - this->energy_llast) < this->of_tole)
        energyConv = true;

    if (this->of_conv == "energy" && energyConv)
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
// Get dL/dphi = dL/drho * drho/dphi = (dE/drho - mu) * 2 * ptempPhi and store it in rdLdphi
// 
void ESolver_OF::calV(double *ptempPhi, double *rdLdphi)
{

    double *dEdtempPhi = new double[this->nrxx];
    double **tempRho = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        tempRho[is] = new double[this->nrxx];
        if (is == this->tnSpinFlag)
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                tempRho[is][ir] = ptempPhi[ir] * ptempPhi[ir];
            }
        }
        else
        {
            for (int ir = 0; ir < this->nrxx; ++ir)
            {
                tempRho[is][ir] = this->ppsi(0, is, ir) * this->ppsi(0, is, ir);
            }   
        }
    }

    GlobalC::pot.vr = GlobalC::pot.v_of_rho(tempRho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    this->tf.get_potential(tempRho);
    for (int i = 0; i < this->nrxx; ++i)
    {
        dEdtempPhi[i] = (GlobalC::pot.vr_eff(this->tnSpinFlag,i) + this->tf.potential[this->tnSpinFlag][i]) * 2. * ptempPhi[i];
    }
    double tempMu = this->cal_mu(ptempPhi, dEdtempPhi, this->nelec[this->tnSpinFlag]);
    for (int i = 0; i < this->nrxx; ++i)
    {
        rdLdphi[i] = dEdtempPhi[i] - 2. * tempMu * ptempPhi[i];
    }
    delete[] dEdtempPhi;
    for (int is = 0; is < GlobalV::NSPIN; ++is) delete[] tempRho[is];
    delete[] tempRho;
} 

//
// Calculate dE/dTheta
// dE/dTheta = <dE/dtempPhi|dtempPhi/dTheta>
//           = <dE/dtempPhi|-phi*sin(theta)+d*cos(theta)>
//
void ESolver_OF::caldEdtheta(double **ptempPhi, double **ptempRho, double *ptheta, double *rdEdtheta)
{
    // double *pdEdtempPhi = new double[this->nrxx];
    double *pdPhidTheta = new double[this->nrxx];

    GlobalC::pot.vr = GlobalC::pot.v_of_rho(ptempRho, GlobalC::CHR.rho_core);
    GlobalC::pot.set_vr_eff();

    this->tf.get_potential(ptempRho);
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nrxx; ++ir)
        {
            this->pdEdphi[is][ir] = (GlobalC::pot.vr_eff(is,ir) + tf.potential[is][ir]) * 2 * ptempPhi[is][ir];

            pdPhidTheta[ir] = - this->ppsi(0, is, ir) * sin(ptheta[is]) + this->pdirect[is][ir] * cos(ptheta[is]);
        }
        rdEdtheta[is] = this->inner_product(this->pdEdphi[is], pdPhidTheta, this->nrxx, this->dV);
    }
    // delete[] pdEdtempPhi;
    delete[] pdPhidTheta;
}

// 
// Calculate chemical potential mu.
// mu = <dE/dphi|phi> / 2nelec.
// 
double ESolver_OF::cal_mu(double *pphi, double *pdEdphi, double nelec)
{
    double mu = this->inner_product(pphi, pdEdphi, this->nrxx, this->dV);
    mu = mu / (2.0*nelec);
    return mu;
}

// 
// print nessecary information
// 
void ESolver_OF::printInfo()
{
    double minDen = GlobalC::CHR.rho[0][0];
    double maxDen = GlobalC::CHR.rho[0][0];
    double minPot = this->pdLdphi[0][0];
    double maxPot = this->pdLdphi[0][0];
    for (int i = 0; i < this->nrxx; ++i)
    {
        if (GlobalC::CHR.rho[0][i] < minDen) minDen = GlobalC::CHR.rho[0][i];
        if (GlobalC::CHR.rho[0][i] > maxDen) maxDen = GlobalC::CHR.rho[0][i];
        if (this->pdLdphi[0][i] < minPot) minPot = this->pdLdphi[0][i];
        if (this->pdLdphi[0][i] > maxPot) maxPot = this->pdLdphi[0][i];
    }
    if (this->iter == 0) cout << "Iter        Etot(Ha)          Theta       PotNorm        min/max(den)          min/max(dL/dPhi)" << endl;
    else cout << setw(6) << this->iter 
    << setw(22) << setiosflags(ios::scientific) << setprecision(12) << this->energy_current/2. 
    << setw(12) << setprecision(3) << this->theta[0]
    << setw(12) << this->normdLdphi
    << setw(10) << minDen << "/ " << setw(12) << maxDen
    << setw(10) << minPot << "/ " << setw(10) << maxPot << endl;
}

void ESolver_OF::cal_Force(ModuleBase::matrix& force)
{
    Forces ff;
    ff.init(force);
}
void ESolver_OF::cal_Stress(ModuleBase::matrix& stress)
{
    Stress_PW ss;
    ss.cal_stress(stress);
    this->tf.get_stress(GlobalC::ucell.omega);
    stress -= this->tf.stress;
}
}