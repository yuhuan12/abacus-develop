#include "H_Hartree_pw.h"
#include "diago_cg.h"
#define eb_k 80
#define EDEPS 1.0
#define ediffsol 1e-3

double H_Hartree_pw::hartree_energy=0.0;

void test_print(double* data, int size)
{
    for(int i=0; i<size; i++)
    {
        cout<<data[i]<<endl;
    }
    return;
}

void test_print(complex<double>* data, int size)
{
    for(int i=0; i<size; i++)
    {
        cout<<data[i].real()<<" "<<data[i].imag()<<endl;
    }
    return;
}

void test_tot_rho(complex<double> *rho, int len, double omega)
{
    double res=0;
    for(int i=0; i<len; i++)
    {
        res+=rho[i].real();
    }
    cout<<"~~~~~~~ test rho ~~~~~~~~"<<endl;
    cout<<res * omega / len<<endl;
    cout<<endl;
    return;
}

void test_tot_rho_G(complex<double> *rho_G, PW_Basis &pwb, double omega)
{
    complex<double> *rho_R = new complex<double>[pwb.nrxx];
    for(int i=0; i<pwb.nrxx; i++)
    {
        rho_R[i].real(0);
        rho_R[i].imag(0);
    }
    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
        rho_R[pwb.ig2fftc[ig]] = rho_G[ig];
    }
    pwb.FFT_chg.FFT3D(rho_R, 1);
    test_tot_rho(rho_R, pwb.nrxx, omega);

}

void H_Hartree_pw::test_res(const UnitCell &ucell,
    PW_Basis &pwb,
    const complex<double>* tot_N,
    complex<double> *phi,
    double* d_eps)
{
    complex<double> *Ax = new complex<double>[pwb.ngmc];

    complex<double> *gradphi_x = new complex<double>[pwb.ngmc];
    complex<double> *gradphi_y = new complex<double>[pwb.ngmc];
    complex<double> *gradphi_z = new complex<double>[pwb.ngmc];

    complex<double> *phi_work = new complex<double>[pwb.ngmc];

    // double *d_eps = new double[pwb.nrxx*2];
    // for(int i=0; i<pwb.nrxx*2; i++)
    // {
    //     d_eps[i] = 1.0;
    // }
    
    Leps(
        ucell,
        pwb,
        phi,
        d_eps,
        gradphi_x,
        gradphi_y,
        gradphi_z,
        phi_work,
        Ax
    );

    // cout <<"~~~~~~~~~~ test ~~~~~~~~~~~~~"<<endl;

    // cout<<"Ax:"<<endl;
    // test_print(Ax, 10);

    // cout<<"b: "<<endl;
    // test_print((complex<double>*)tot_N, 10);

    // cout<<"============= result phi: ============="<<endl;
    // test_print(phi, 10);
}

//--------------------------------------------------------------------
// Transform charge density to hartree potential.
//--------------------------------------------------------------------
ModuleBase::matrix H_Hartree_pw::v_hartree(
	const UnitCell &cell, 
	PW_Basis &pwb, 
	const int &nspin,
	const double*const*const rho)
{
    ModuleBase::TITLE("H_Hartree_pw","v_hartree");
    ModuleBase::timer::tick("H_Hartree_pw","v_hartree");

    //  Hartree potential VH(r) from n(r)
    std::vector<std::complex<double>> Porter(pwb.nrxx);
    const int nspin0 = (nspin==2) ? 2 : 1;
    for(int is=0; is<nspin0; is++)
        for (int ir=0; ir<pwb.nrxx; ir++) 
            Porter[ir] += std::complex<double>( rho[is][ir], 0.0 );
    //=============================
    //  bring rho (aux) to G space
    //=============================
    pwb.FFT_chg.FFT3D(Porter.data(), -1);

    //double charge;
    //if (pwb.gstart == 1)
    //    charge = cell.omega * Porter[pwb.ig2fftc[0]].real();
    //OUT(GlobalV::ofs_running, "v_h charge", charge);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================

	double ehart = 0.0;

    std::vector<std::complex<double>> vh_g(pwb.ngmc);

    // cout<<"Porter_g in v_hartree"<<endl;
    // cout<<"pwb.gstart:"<<pwb.gstart<<endl;
    for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * pwb.gg [ig]);

            ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
            vh_g[ig] = fac * Porter[j];
            // if(ig<10)
            // {
            //     cout<<Porter[j].real()<<" "<<Porter[j].imag()<<endl;
            // }
        }
    }

    Parallel_Reduce::reduce_double_pool( ehart );
    ehart *= 0.5 * cell.omega;
    //std::cout << " ehart=" << ehart << std::endl;
    H_Hartree_pw::hartree_energy = ehart;

    std::fill( Porter.begin(), Porter.end(), std::complex<double>(0.0,0.0) );
    for (int ig = 0;ig < pwb.ngmc;ig++)
        Porter[pwb.ig2fftc[ig]] = vh_g[ig];
    //==========================================
    //transform hartree potential to real space
    //==========================================
    pwb.FFT_chg.FFT3D(Porter.data(), 1);

    //==========================================
    //Add hartree potential to the xc potential
    //==========================================
    ModuleBase::matrix v(nspin, pwb.nrxx);
    if(nspin==4)
    {
        for (int ir = 0;ir < pwb.nrxx;ir++)
            v(0, ir) = Porter[ir].real();
    }
    else
    {
        for (int is = 0;is < nspin;is++)
            for (int ir = 0;ir < pwb.nrxx;ir++)
                v(is, ir) = Porter[ir].real();
    }

//-----------------------------------------------------------
// we need to add this out_potential funciton back 
// in near future, 2021-02-25
//-----------------------------------------------------------
	//-------------------------------------------
	// output the Hartree potential into a file.
	//-------------------------------------------
/*
	if(out_potential==-2)
	{
		std::cout << " output VH" << std::endl;
		int is = 0;
		int iter = 0;
		int precision = 3;
		std::string fn = "VH.dat";
		std::stringstream ss;
		ss << GlobalV::global_out_dir << fn;
		ModuleBase::matrix v;
		v.create(1,pwb.nrxx);
		for(int ir=0; ir<pwb.nrxx; ++ir)
		{
			v(0,ir) = Porter[ir].real();
		}
		this->write_potential( is, iter, ss.str(), v, precision, 1 );
	}
*/

    ModuleBase::timer::tick("H_Hartree_pw","v_hartree");
    return v;
} // end subroutine v_h

void H_Hartree_pw::gauss_charge(const UnitCell &cell, PW_Basis &pwb, complex<double> *N, const int flag)
{
	const double delta_grd = pow(cell.omega / pwb.ngmc, 1.0 / 3) ;    //unit bohr??
    const double sigma_nc_k = 1.6 * delta_grd ;
	
	// complex<double> *N = new complex<double>[pwb.ngmc] ;

    pwb.setup_structure_factor();   //call strucFac(ntype,ngmc) 

    for (int it=0; it<cell.ntype; it++)
    {
	    double RCS = cell.atoms[it].r[cell.atoms[it].mesh-1] ;
	    double sigma_rc_k = RCS / 2.5 ;

        for (int ia=0; ia<cell.atoms[it].na; ia++)
		{
            for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)  // ngmc : num. of G vectors within ggchg in each proc.
			{
                //G^2
                double gg = pwb.get_NormG_cartesian(ig);
		        gg = gg * cell.tpiba2; 

                // gauss ionic charge
				if(flag == 1)    
                {
                    //need Z(proton num.) 
                    N[ig].real(N[ig].real() - get_Z(cell.atoms[it].psd) * pwb.strucFac(it, ig).real() * exp(-0.5 *gg*(sigma_nc_k * sigma_nc_k)));
                    N[ig].imag(N[ig].imag() - get_Z(cell.atoms[it].psd) * pwb.strucFac(it, ig).imag() * exp(-0.5 *gg*(sigma_nc_k * sigma_nc_k)));
				}
				if(flag == 3)    // pseudo core
                {
                    N[ig].real(N[ig].real()+(get_Z(cell.atoms[it].psd)-cell.atoms[it].zv) * pwb.strucFac(it, ig).real() * exp(-0.5*gg*(sigma_rc_k*sigma_rc_k))); 
                    N[ig].imag(N[ig].imag()+(get_Z(cell.atoms[it].psd)-cell.atoms[it].zv) * pwb.strucFac(it, ig).imag() * exp(-0.5*gg*(sigma_rc_k*sigma_rc_k))); 
				}
            }
            // cout<<"test getZ"<<endl;
            // cout<<get_Z(cell.atoms[it].psd)<<endl;
            // cout<<cell.atoms[it].zv<<endl;
 		}
    }
}

int H_Hartree_pw::get_Z(string str)
{
    if(str == "H") 
    {
        return 1;
    }
    else if(str == "O") 
    {
        return 8;
    }
    else return 1;
}

//

ModuleBase::matrix H_Hartree_pw::v_correction(const UnitCell &cell,
                                PW_Basis &pwb,
                                const int &nspin, 
                                const double*const*const rho)
{				
    ModuleBase::TITLE("H_Hartree_pw","v_correction");
    ModuleBase::timer::tick("H_Hartree_pw","v_correction");	
    // Porter in real space	
    // Get sum on nspin axis .	
    complex<double> *Porter = new complex<double>[pwb.nrxx];
    const int nspin0 = (nspin==2) ? 2 : 1;
    for(int is=0; is<nspin0; is++)
        for (int ir=0; ir<pwb.nrxx; ir++) 
            Porter[ir] += std::complex<double>( rho[is][ir], 0.0 );
    
    cout<< "*********** abacus origin: ***********"<<endl;
    test_tot_rho(Porter, pwb.nrxx, cell.omega);
    
    // into G space (nrxx)
    pwb.FFT_chg.FFT3D(Porter, -1);

    // Build Porter_g in G space. (dim = ngmc)
    complex<double> *Porter_g = new complex<double>[pwb.ngmc];
    ModuleBase::GlobalFunc::ZEROS( Porter_g, pwb.ngmc);
    for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            Porter_g[ig] = Porter[j];
        }
    }
    // delete Porter in real space
    delete [] Porter;

    cout<< "*********** Porter G (origin): ***********"<<endl;
    test_tot_rho_G(Porter_g, pwb, cell.omega);

    // From now, Porter_g is in G space (dim=ngmc)
    // TOTN(n_val+N_gauss)
    complex<double> *N = new complex<double>[pwb.ngmc];
    complex<double> *TOTN = new complex<double>[pwb.ngmc];
    complex<double> *PS_TOTN = new complex<double>[pwb.ngmc];

    // N in G space
    gauss_charge(cell, pwb, N, 1); 

    cout<< "*********** N (gauss param = 1): ***********"<<endl;
    test_tot_rho_G(N, pwb, cell.omega);
    
    // TOTN in G space
    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
		TOTN[ig] = N[ig] + Porter_g[ig] ;
    }

    cout<< "*********** TOTN (N + origin): ***********"<<endl;
    test_tot_rho_G(TOTN, pwb, cell.omega);

    // PS_TOTN(n_val+pseudo_core)
    gauss_charge(cell, pwb, N, 3);

    cout<< "*********** N (gauss param = 3): ***********"<<endl;
    test_tot_rho_G(N, pwb, cell.omega);

    // PS_TOTN in G space (ngmc)
    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
		PS_TOTN[ig] = N[ig] + Porter_g[ig] ;
    }

    cout<< "*********** PS_TOT_N before FFT: ***********"<<endl;
    test_tot_rho_G(PS_TOTN, pwb, cell.omega);

    // Build a nrxx vector to DO FFT .
    complex<double> *PS_TOTN_real = new complex<double>[pwb.nrxx];
    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
        PS_TOTN_real[pwb.ig2fftc[ig]] = PS_TOTN[ig];
    }

    pwb.FFT_chg.FFT3D(PS_TOTN_real, 1);
    
    cout<< "*********** PS_TOT_N after FFT: ***********"<<endl;
    test_tot_rho(PS_TOTN_real, pwb.nrxx, cell.omega);

    // shapefunction value varies from 0 in the solute to 1 in the solvent
    // epsilon = 1.0 + (eb_k - 1) * shape function

    double nc_k = 0.00037; //ang^-3
    double sigma_k = 0.6 ;

    // build epsilon in real space (nrxx)
    double *epsilon = new double[pwb.nrxx];
    double *shapefunc = new double[pwb.nrxx];
    for(int i=0; i<pwb.nrxx; i++)
    {
	    shapefunc[i] = erfc((log(max(PS_TOTN_real[i].real(), 1e-10) / nc_k)) / sqrt(2.0) / sigma_k ) / 2;
		epsilon[i] = 1 + (eb_k - 1) * shapefunc[i];
	}

    // cout<<"build shapefunc"<<endl;
    // test_print(shapefunc, 15);

    // cout<<"build epsilon"<<endl;
    // test_print(epsilon+50000, 10000);

    delete [] PS_TOTN_real;

    // fill in epsilon with ones (Used to test)
    // for(int i=0; i<pwb.nrxx; i++)
    // {
    //     epsilon[i] = 1.0;
    // }

    // test_print(epsilon+50000, 100);

    complex<double> *Sol_phi = new complex<double>[pwb.ngmc];
    int ncgsol = 0;

    // call CG method
    cout << "call minimize" << endl;
    cout << "ngmc = " << pwb.ngmc << endl;
    cout << "nrxx = " << pwb.nrxx << endl;
    minimize(cell, pwb, epsilon, TOTN, Sol_phi, ncgsol);

    // fake v
    ModuleBase::matrix v(nspin, pwb.nrxx);

    delete [] Porter_g;
    delete [] Sol_phi;
    delete [] N;
    delete [] PS_TOTN;
    delete [] TOTN;
    // delete [] r_work;
    delete [] epsilon;
    delete [] shapefunc;

    ModuleBase::timer::tick("H_Hartree_pw","v_correction");	
    return v;

}


// cast complex to real
void H_Hartree_pw::cast_C2R(complex<double> *src, double* dst, int dim)
{
    int k = 0;
    for(int i=0; i<dim; i++)
    {
        dst[k] = src[i].real();
        dst[k+1] = src[i].imag();
        k = k + 2;
    }
    return ;
}


// Routine for solving poissons eqn using conjugate gradient (CG Method) ..
void H_Hartree_pw::minimize(
    const UnitCell &ucell,
    PW_Basis &pwb,
    double *d_eps, // dim=nrxx
    const complex<double>* tot_N,
    complex<double> *phi,
    int &ncgsol // output
)
{
    // parameters of CG method
    double alpha = 0;
    double beta = 0;
    // r * r'
    double rinvLr = 0;
    // r * r
    double r2 = 0;
    // precond loop parameter
    int i = 0;
    // malloc vectors in G space 
    complex<double> *resid = new complex<double>[pwb.ngmc];
    complex<double> *z = new complex<double>[pwb.ngmc];
    complex<double> *lp = new complex<double>[pwb.ngmc];
    complex<double> *gsqu = new complex<double>[pwb.ngmc];
    complex<double> *d = new complex<double>[pwb.ngmc];

    complex<double> *gradphi_x = new complex<double>[pwb.ngmc];
    complex<double> *gradphi_y = new complex<double>[pwb.ngmc];
    complex<double> *gradphi_z = new complex<double>[pwb.ngmc];

    complex<double> *phi_work = new complex<double>[pwb.ngmc];

    int count = 0;
    double gg = 0;

    // calculate precondition vector GSQU (In G space, ngmc)
    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
        gg = pwb.get_NormG_cartesian(ig);
        gsqu[ig].real(1.0 / (gg * ucell.tpiba2 )); // without kappa_2
        gsqu[ig].imag(0);
    }

    // cout<<"precond"<<endl;
    // test_print(gsqu, 10);

    // init guess for phi
    for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        // sol_phi = ?
        phi[ig].real(tot_N[ig].real() * gsqu[ig].real());
        phi[ig].imag(tot_N[ig].imag() * gsqu[ig].real()); 

        // phi[ig].real(0.0);
        // phi[ig].imag(0.0);
    }
    
    // call leps to calculate div ( epsilon * grad ) phi - kappa^2 * phi

    cout<<"before CG"<<endl;
    test_res(ucell, pwb, tot_N, phi, d_eps);

    // cout <<"****phi before"<<endl;
    // test_print(phi, 20);
    Leps(
        ucell,
        pwb,
        phi,
        d_eps,
        gradphi_x,
        gradphi_y,
        gradphi_z,
        phi_work,
        lp
    );

    // int c;
    // cin >>c;

    // cout <<"****phi after"<<endl;
    // test_print(phi, 10);

    cout << "------- lp ---------"<<endl;
    test_print(lp, 10);
    // calculate Lp: pwb.ngmc, complex-double

    // the residue
    // r = A*phi + (chtot + N)
    // pwb.ngmc = ?
    for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        // resid[ig].real(lp[ig].real() + tot_N[ig].real());
        // resid[ig].imag(lp[ig].imag() + tot_N[ig].imag());
        resid[ig] = lp[ig] + tot_N[ig];
    }

    // cout << "-------resid now-----------"<<endl;
    // test_print(resid, 10);

    cout << "r2:" <<endl;
    cout<< Diago_CG::ddot_real(pwb.ngmc, resid, resid)<<endl;
    // precondition of the residue, z = invLr
    for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        z[ig].real(gsqu[ig].real() * resid[ig].real());
        z[ig].imag(gsqu[ig].real() * resid[ig].imag());
    }
    // calculate r*r'
    // where ddot_real()
    // Diago_CG::ddot_real( npw, respsi, respsi );
    rinvLr = Diago_CG::ddot_real(pwb.ngmc, resid, z);
    r2 = Diago_CG::ddot_real(pwb.ngmc, resid, resid);

    double r20 = r2;

    // copy
    for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        // d[ig].real(z[ig].real());
        // d[ig].imag(z[ig].imag());
        d[ig] = z[ig];
    }

    // CG Loop
    cout << "r2: " << r2 << endl;
    cout << "sqrt(r2): " << sqrt(r2) << endl;
    cout << "sqrt(rinvLr): " << sqrt(rinvLr) << endl;
    while(count < 20000 && sqrt(r2) > ediffsol && sqrt(rinvLr) > 1e-10)
    {
        cout << "=====================  CG Loop: " << count <<"===================" <<endl;
        cout << "sqrt(r2): " << sqrt(r2) << endl;
        cout << "sqrt(rinvLr): " << sqrt(rinvLr) << endl;

        if(sqrt(r2)>1e6)
        {
            cout << "CG ERROR!!!" <<endl;
            break;
        }
        
        Leps(
            ucell,
            pwb,
            d,
            d_eps,
            gradphi_x,
            gradphi_y,
            gradphi_z,
            phi_work,
            lp
        );

        // cout <<"lp after leps"<<endl;
        // test_print(lp, 10);
        // calculate alpha
        alpha = -rinvLr / Diago_CG::ddot_real(pwb.ngmc, d, lp);
        cout<<"alpha: "<<alpha<<endl;
        // update phi
        for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            // phi[ig].real() += alpha * d[ig].real();
            // phi[ig].imag() += alpha * d[ig].imag();
            phi[ig] += alpha * d[ig];
        }
        // update resid

        // cout << "-------resid before update -----------"<<endl;
        // test_print(resid, 10);
        for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            // resid[ig].real() += alpha * lp[ig].real();
            // resid[ig].imag() += alpha * lp[ig].imag();
            resid[ig] += alpha * lp[ig];
        }

        // cout << "-------resid after update -----------"<<endl;
        // test_print(resid, 10);
        // precond one more time..
        for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            // z[ig].real() = gsqu[ig].real() * resid[ig].real();
            // z[ig].imag() = gsqu[ig].real() * resid[ig].imag();
            z[ig] = gsqu[ig] * resid[ig];
        }

        cout << "r2:" <<endl;
        cout<< Diago_CG::ddot_real(pwb.ngmc, resid, resid)<<endl;
        // calculate beta
        beta = 1.0 / rinvLr;
        rinvLr = Diago_CG::ddot_real(pwb.ngmc, resid, z);
        beta *= rinvLr;
        // update d
        for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
        {
            // d[ig].real() = beta * d[ig].real() + z[ig].real();
            // d[ig].imag() = beta * d[ig].imag() + z[ig].imag();
            d[ig] = beta * d[ig] + z[ig];
        }
        r2 = 0;
        r2 = Diago_CG::ddot_real(pwb.ngmc, resid, resid);

        // update counter
        count++;
    }// end CG loop

    cout<<"after CG"<<endl;
    test_res(ucell, pwb, tot_N, phi, d_eps);
    // output: num of cg loop
    ncgsol = count;
    
    // TODO:
    for(int ig = pwb.gstart; ig < pwb.ngmc; ig++)
    {
        // phi[ig].real() = phi[ig].real() * EDEPS / ucell.omega;
        // phi[ig].imag() = phi[ig].imag() * EDEPS / ucell.omega;
        phi[ig] = phi[ig] * EDEPS / ucell.omega;
    }
    // free
    delete [] resid;
    delete [] z;
    delete [] lp;
    delete [] gsqu;
    delete [] d;
    delete [] gradphi_x;
    delete [] gradphi_y;
    delete [] gradphi_z;
    delete [] phi_work;
}


void H_Hartree_pw::Leps(
    const UnitCell &ucell,
    PW_Basis &pwb,
    complex<double> *phi,
    double *epsilon, // epsilon from shapefunc, dim=nrxx
    complex<double> *gradphi_x, // dim=ngmc
    complex<double> *gradphi_y,
    complex<double> *gradphi_z,
    complex<double> *phi_work,
    complex<double> *lp // output
)
{
    for(int i=pwb.gstart; i<pwb.ngmc; i++)
    {
        gradphi_x[i].real(phi[i].imag() * pwb.gcar[i].x * ModuleBase::TWO_PI);
        gradphi_y[i].real(phi[i].imag() * pwb.gcar[i].y * ModuleBase::TWO_PI);
        gradphi_z[i].real(phi[i].imag() * pwb.gcar[i].z * ModuleBase::TWO_PI);

        gradphi_x[i].imag(-phi[i].real() * pwb.gcar[i].x * ModuleBase::TWO_PI);
        gradphi_y[i].imag(-phi[i].real() * pwb.gcar[i].y * ModuleBase::TWO_PI);
        gradphi_z[i].imag(-phi[i].real() * pwb.gcar[i].z * ModuleBase::TWO_PI);

        // gradphi_x[i].real(phi[i].imag() * pwb.gcar[i].x * 1);
        // gradphi_y[i].real(phi[i].imag() * pwb.gcar[i].y * 1);
        // gradphi_z[i].real(phi[i].imag() * pwb.gcar[i].z * 1);

        // gradphi_x[i].imag(-phi[i].real() * pwb.gcar[i].x * 1);
        // gradphi_y[i].imag(-phi[i].real() * pwb.gcar[i].y * 1);
        // gradphi_z[i].imag(-phi[i].real() * pwb.gcar[i].z * 1);
    }
    // build real space vectors todo FFT
    complex<double> *gradphi_x_real = new complex<double>[pwb.nrxx];
    complex<double> *gradphi_y_real = new complex<double>[pwb.nrxx];
    complex<double> *gradphi_z_real = new complex<double>[pwb.nrxx];

    complex<double> *phi_real = new complex<double>[pwb.nrxx];
    complex<double> *phi_work_real = new complex<double>[pwb.nrxx];

    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
        gradphi_x_real[pwb.ig2fftc[ig]] = gradphi_x[ig];
        gradphi_y_real[pwb.ig2fftc[ig]] = gradphi_y[ig];
        gradphi_z_real[pwb.ig2fftc[ig]] = gradphi_z[ig];

        phi_real[pwb.ig2fftc[ig]] = phi[ig];
    }

    // fft into real space
    pwb.FFT_chg.FFT3D(gradphi_x_real, 1);
    pwb.FFT_chg.FFT3D(gradphi_y_real, 1);
    pwb.FFT_chg.FFT3D(gradphi_z_real, 1);

    pwb.FFT_chg.FFT3D(phi_real, 1);

    // cout<<"epsilon"<<endl;
    // test_print(epsilon, 10);

    // multiple grad x, y, z with epsilon in real space(nrxx) 
    for(int j=0; j<pwb.nrxx; j++)
    {
        // gradphi_x_real[j].real(gradphi_x_real[j].real() * epsilon[2*j]);
        // gradphi_y_real[j].real(gradphi_y_real[j].real() * epsilon[2*j]);
        // gradphi_z_real[j].real(gradphi_z_real[j].real() * epsilon[2*j]);

        // gradphi_x_real[j].imag(gradphi_x_real[j].imag() * epsilon[2*j+1]);
        // gradphi_y_real[j].imag(gradphi_y_real[j].imag() * epsilon[2*j+1]);
        // gradphi_z_real[j].imag(gradphi_z_real[j].imag() * epsilon[2*j+1]);

        // dim of epsilon = pwb.nrxx * 1 ! (complex mode)
        gradphi_x_real[j] *= epsilon[j];
        gradphi_y_real[j] *= epsilon[j];
        gradphi_z_real[j] *= epsilon[j];

        phi_work_real[j] = phi_real[j] * eb_k;
    }

    // FFT inverse
    pwb.FFT_chg.FFT3D(gradphi_x_real, -1);
    pwb.FFT_chg.FFT3D(gradphi_y_real, -1);
    pwb.FFT_chg.FFT3D(gradphi_z_real, -1);

    pwb.FFT_chg.FFT3D(phi_work_real, -1);
    pwb.FFT_chg.FFT3D(phi_real, -1);

    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
        gradphi_x[ig] = gradphi_x_real[pwb.ig2fftc[ig]];
        gradphi_y[ig] = gradphi_y_real[pwb.ig2fftc[ig]];
        gradphi_z[ig] = gradphi_z_real[pwb.ig2fftc[ig]];

        phi_work[ig] = phi_work_real[pwb.ig2fftc[ig]];
        phi[ig] = phi_real[pwb.ig2fftc[ig]];
    }

    delete [] gradphi_x_real;
    delete [] gradphi_y_real;
    delete [] gradphi_z_real;

    delete [] phi_work_real;
    delete [] phi_real;

    // div(epsilon*grad phi) in kspace
    // add the kappa^2 contrib, -phi_work
    for(int ig=pwb.gstart; ig<pwb.ngmc; ig++)
    {
        // use 2*pi or not ?
        // lp[ig] = ModuleBase::TWO_PI * 
        //     (gradphi_x[ig]*pwb.gcar[ig].x + gradphi_y[ig]*pwb.gcar[ig].y + gradphi_z[ig]*pwb.gcar[ig].z
        //     - phi_work[ig]);
        
        double tmp_real = gradphi_x[ig].real()*pwb.gcar[ig].x + 
            gradphi_y[ig].real()*pwb.gcar[ig].y + 
            gradphi_z[ig].real()*pwb.gcar[ig].z;

        double tmp_imag = gradphi_x[ig].imag()*pwb.gcar[ig].x + 
            gradphi_y[ig].imag()*pwb.gcar[ig].y + 
            gradphi_z[ig].imag()*pwb.gcar[ig].z;
        
        lp[ig].real(ModuleBase::TWO_PI * tmp_imag);
        lp[ig].imag(-ModuleBase::TWO_PI * tmp_real);

        // use phi_work ?
        // lp[ig] -= phi_work[ig];

    }
}

/*
void H_Hartree_pw::createcavity(
    const UnitCell &ucell,
    PW_Basis &pwb,
    const complex<double>* PS_TOTN  
)
{

    
    ModuleBase::Vector3<double> *nablan ;
    double *nablan_2 = (double*)malloc(pwb.ngmc * sizeof(double));
    double *sqrt_nablan_2 = (double*)malloc(pwb.ngmc * sizeof(double));
    double *lapn = (double*)malloc(pwb.ngmc * sizeof(double));

    // nabla n
    GGA_PW::grad_rho( PS_TOTN , nablan );

    //  |\nabla n |^2 = nablan_2 
    for (int ir=0;ir<pwb.nrxx;ir++)
    {
        nablan_2[ir] = pow(nablan[ir].x,2) +pow(nablan[ir].y,2)+pow(nablan[ir].z,2) ;
    }
    
    // Laplacian of n
    lapl_rho( PS_TOTN , lapn) ;

    //-------------------------------------------------------------
    //add -Lap(n)/|\nabla n| to vwork and copy \sqrt(|\nabla n|^2) 
    //to sqrt_nablan_2
    //-------------------------------------------------------------
    double *vwork = (double*)malloc(pwb.ngmc * sizeof(double));

    for (int ir=0;ir<pwb.nrxx;ir++)
    {
        tmp = sqrt(nablan_2[ir]);
        vwork(ir) = vwork(ir) - (lapn[ir])/tmp ;
        sqrt_nablan_2[ir] = tmp ;
    }

    //-------------------------------------------------------------
    //term1 = gamma*A / n, where
    //gamma * A = exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )
    //-------------------------------------------------------------
    shape_gradn(PS_TOTN , term1) ;

    //-------------------------------------------------------------
    //quantum surface area, integral of (gamma*A / n) * |\nabla n| 
    //=term1 * sqrt_nablan_2
    //-------------------------------------------------------------
    double qs = 0;
    double Acav = 0;

    for (int ir=0;ir<pwb.nrxx;ir++)
    {
        qs = qs + (term1[ir]) * (sqrt_nablan_2[ir]) ;

        //   1/ |nabla n|
        sqrt_nablan_2[ir] = 1/sqrt_nablan_2[ir] ;

    }
    //-------------------------------------------------------------
    //cavitation energy
    //-------------------------------------------------------------
    Acav = tau * qs ;

    //  packs the real array into a complex one
    complex<double> *inv_gn = new complex<double>[pwb.nrxx];
    cast_R2C(sqrt_nablan_2, inv_gn, pwb.nrxx);

    //  to G space
    pwb.FFT_chg.FFT3D(inv_gn.data(), -1);
    complex<double> *inv_gn_g = new complex<double>[pwb.ngmc];

     for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) 
        {
            inv_gn_g[ig] = inv_gn[j] ;
        }
    }

    // \nabla(1 / |\nabla n|), ggn in real space
    ModuleBase::Vector3<double> *ggn;
    GGA_PW::grad_rho( inv_gn_g , ggn );

    //-------------------------------------------------------------
    //add -(\nabla n . \nabla(1/ |\nabla n|)) to Vcav in real space
    //and multiply by term1 = gamma*A/n in real space
    //-------------------------------------------------------------
    for(int ir=0;ir<pwb.nrxx;ir++)
    {
        tmp = (nablan[ir].x*ggn[ir].x + nablan[ir].y*ggn[ir].y + nablan[ir].z*ggn[ir].z) * term1[ir] ;
        vwork[ir] = vwork[ir] - tmp;
    }

    complex<double> *Vcav = new complex<double>[pwb.nrxx];
    cast_R2C(vwork, Vcav, pwb.nrxx) ;

    //  to G space
    pwb.FFT_chg.FFT3D(Vcav.data(), -1);
    complex<double> *Vcav_g = new complex<double>[pwb.ngmc];

     for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) 
        {
            Vcav_g[ig] = Vcav[j] ;
        }
    }
    Vcav_g = Vcav_g *tau ;



}


void H_Hartree_pw::lapl_rho( const std::complex<double> *rhog, double *lapn )
{
    std::complex<double> *gdrtmpg = new std::complex<double>[GlobalC::pw.ngmc];
    ModuleBase::GlobalFunc::ZEROS(gdrtmpg, GlobalC::pw.ngmc);

    std::complex<double> *Porter = GlobalC::UFFT.porter;

    // the formula is : rho(r)^prime = \int iG * rho(G)e^{iGr} dG
    for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
        gdrtmpg[ig] = rhog[ig];

    // calculate the charge density gradient in reciprocal space.
    ModuleBase::GlobalFunc::ZEROS(Porter, GlobalC::pw.nrxx);
    for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
        Porter[ GlobalC::pw.ig2fftc[ig] ] = gdrtmpg[ig]* pow(std::complex<double>(GlobalC::pw.get_G_cartesian_projection(ig, 0), 0.0),2);
    // bring the gdr from G --> R
    GlobalC::pw.FFT_chg.FFT3D(Porter, 1);
    // remember to multily 2pi/a0, which belongs to G vectors.
    for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
        lapn[ir] =lapn[ir]-Porter[ir].real() * pow(GlobalC::ucell.tpiba，2）;

    // calculate the charge density gradient in reciprocal space.
    ModuleBase::GlobalFunc::ZEROS(Porter, GlobalC::pw.nrxx);
    for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
        Porter[GlobalC::pw.ig2fftc[ig]] = gdrtmpg[ig] * pow(std::complex<double>(GlobalC::pw.get_G_cartesian_projection(ig, 1), 0.0),2);
    // bring the gdr from G --> R
    GlobalC::pw.FFT_chg.FFT3D(Porter, 1);
    // remember to multily 2pi/a0, which belongs to G vectors.
    for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
        lapn[ir] = lapn[ir]-Porter[ir].real() * pow(GlobalC::ucell.tpiba,2);

    // calculate the charge density gradient in reciprocal space.
    ModuleBase::GlobalFunc::ZEROS(Porter, GlobalC::pw.nrxx);
    for(int ig=0; ig<GlobalC::pw.ngmc; ig++)
        Porter[GlobalC::pw.ig2fftc[ig]] = gdrtmpg[ig] * pow(std::complex<double>(GlobalC::pw.get_G_cartesian_projection(ig, 2), 0.0),2);
    // bring the gdr from G --> R
    GlobalC::pw.FFT_chg.FFT3D(Porter, 1);
    // remember to multily 2pi/a0, which belongs to G vectors.
    for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
        lapn[ir] = lapn[ir]-Porter[ir].real() * pow(GlobalC::ucell.tpiba,2);

    delete[] gdrtmpg;
    return;
}



//calculates first derivative of the shape function  wrt CHTOT 
//in realspace
// exp(-(log(n/n_c))^2 /(2 sigma^2)) /(sigma * sqrt(2*pi) )/n
void H_Hartree_pw::shape_gradn( const complex<double> *PS_TOTN ,
                                PW_Basis &pw,
                                double *eprime )
{
     // Build a nrxx vector to DO FFT .
    complex<double> *PS_TOTN_real = new complex<double>[pwb.nrxx];

    for(int ig=0; ig<pwb.ngmc; ig++)
    {
        PS_TOTN_real[pwb.ig2fftc[ig]] = PS_TOTN[ig];
    }

    // to real space
    pwb.FFT_chg.FFT3D(PS_TOTN_real, 1);

    double epr_c = 1.0/sqrt(ModuleBase::TWO_PI)/sigma_k ;
    double epr_z = 0 ;
    for(int ir=0; ir<pw.nrxx; ir++)
    {
        epr_z = log(PS_TOTN_real[ir].real()/nc_k)/sqrt(2)/sigma_k ;
        eprime[ir] = epr_c*(exp(-pow(epr_z,2))/PS_TOTN_real[ir].real() ;
    }
    
}
*/
