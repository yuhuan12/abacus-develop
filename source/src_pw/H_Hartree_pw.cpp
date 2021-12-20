#include "H_Hartree_pw.h"
#include "diago_cg.h"
#define eb_k 80.0
#define EDEPS 1.0
#define ediffsol 0.04

double H_Hartree_pw::hartree_energy=0.0;

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
    for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI / (cell.tpiba2 * pwb.gg [ig]);

            ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
            vh_g[ig] = fac * Porter[j];
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
	const double delta_grd = pow(cell.omega / pwb.ngmc, 1.0 / 3) ;    //(GRIDC%NPLWV) = ngmc?
    const double sigma_nc_k = 1.6 * delta_grd ;
	
	// complex<double> *N = new complex<double>[pwb.ngmc] ;

    pwb.setup_structure_factor();   //call strucFac(ntype,ngmc) 

    for (int it=0; it<cell.ntype; it++)
    {
	    double RCS = cell.atoms[it].r[cell.atoms[it].mesh-1] ;
	    double sigma_rc_k = RCS / 2.5 ;

        for (int ia=0; ia<cell.atoms[it].na; ia++)
		{
            for(int ig=0; ig<pwb.ngmc; ig++)  // ngmc : num. of G vectors within ggchg in each proc.
			{
                //G^2
                double gg = pwb.get_NormG_cartesian(ig);
		        gg = gg * ModuleBase::TWO_PI * ModuleBase::TWO_PI; 

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
 		}
    }
}

int H_Hartree_pw::get_Z(string str)
{
    if(str == "H") return 1;
    else if(str == "O") return 8;
    else return 1;
}

//

void H_Hartree_pw::V_correction(const UnitCell &cell,
                                PW_Basis &pwb, 
                                double **rho)
// Porter(nrxx): rho to G space;  call pwb.ig2fftc to get Poter_g(ngmc)
{							
    complex<double> *Porter_g = new complex<double>(pwb.ngmc);
    for (int ig = pwb.gstart; ig<pwb.ngmc; ig++)
    {
        const int j = pwb.ig2fftc[ig];
        if(pwb.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            Porter_g[ig].real(rho[0][j]);
            Porter_g[ig].imag(0.0);
        }
    }
    // TOTN(n_val+N_gauss)
    complex<double> *N = new complex<double>[pwb.ngmc];
    complex<double> *TOTN = new complex<double>[pwb.ngmc];
    complex<double> *PS_TOTN = new complex<double>[pwb.ngmc];
    gauss_charge(cell, pwb, N, 1); 

    for(int ig=0; ig<pwb.ngmc; ig++)
    {
		TOTN[ig] = N[ig] + Porter_g[ig] ;
    }

    // PS_TOTN(n_val+pseudo_core)
    gauss_charge(cell, pwb, N, 3); 

    for(int ig=0; ig<pwb.ngmc; ig++)
    {
		PS_TOTN[ig] = N[ig] + Porter_g[ig] ;
    }

    double *r_work = (double*)malloc(pwb.ngmc * sizeof(double));
    cast_C2R(PS_TOTN, r_work, pwb.ngmc);

    // shapefunction value varies from 0 in the solute to 1 in the solvent
    // epsilon = 1.0 + (eb_k - 1) * shape function

    double nc_k = 0.0025 ; //ang^-3
    double sigma_k = 0.6 ;

    double *epsilon = (double*)malloc(pwb.ngmc * sizeof(double));
    double *shapefunc = (double*)malloc(pwb.ngmc * sizeof(double));
    for(int i=0; i<pwb.ngmc; i++)
    {
	    shapefunc[i] = erfc((log(r_work[i]) / nc_k) / sqrt(2.0) / sigma_k ) / 2;
		epsilon[i] = 1 + (eb_k - 1) * shapefunc[i];
	}

    complex<double> *Sol_phi = new complex<double>[pwb.ngmc];
    int ncgsol = 0;
    minimize(cell, pwb, epsilon, TOTN, Sol_phi, ncgsol);

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

void H_Hartree_pw::Leps(
    const UnitCell &ucell,
    PW_Basis &pwb,
    complex<double> *phi,
    double *epsilon, // epsilon from shapefunc
    complex<double> *gradphi_x,
    complex<double> *gradphi_y,
    complex<double> *gradphi_z,
    complex<double> *phi_work,
    complex<double> *lp // output
)
{
    int n = pwb.ngmc;
    // TODO
    for(int i=0; i<pwb.ngmc; i++)
    {
        gradphi_x[i].real(- phi[i].imag() * pwb.gcar[i].x * ModuleBase::TWO_PI);
        gradphi_y[i].real(- phi[i].imag() * pwb.gcar[i].y * ModuleBase::TWO_PI);
        gradphi_z[i].real(- phi[i].imag() * pwb.gcar[i].z * ModuleBase::TWO_PI);

        gradphi_x[i].imag(phi[i].real() * pwb.gcar[i].x * ModuleBase::TWO_PI);
        gradphi_y[i].imag(phi[i].real() * pwb.gcar[i].y * ModuleBase::TWO_PI);
        gradphi_z[i].imag(phi[i].real() * pwb.gcar[i].z * ModuleBase::TWO_PI);
    }
    // FFT forward
    // todo: FFT1d
    pwb.FFT_chg.FFT3D(gradphi_x, 1);
    pwb.FFT_chg.FFT3D(gradphi_y, 1);
    pwb.FFT_chg.FFT3D(gradphi_z, 1);

    pwb.FFT_chg.FFT3D(phi, 1);

    // calculate
    for(int j=0; j<pwb.ngmc; j++)
    {
        gradphi_x[j] *= epsilon[j];
        gradphi_y[j] *= epsilon[j];
        gradphi_z[j] *= epsilon[j];

        // gradphi_x[j].imag() *= epsilon[j];
        // gradphi_y[j].imag() *= epsilon[j];
        // gradphi_z[j].imag() *= epsilon[j];

        phi_work[j].real(phi[j].real() * eb_k);
        phi_work[j].imag(phi[j].imag() * eb_k);
    }

    // FFT inverse
    pwb.FFT_chg.FFT3D(gradphi_x, -1);
    pwb.FFT_chg.FFT3D(gradphi_y, -1);
    pwb.FFT_chg.FFT3D(gradphi_z, -1);

    pwb.FFT_chg.FFT3D(phi_work, -1);
    pwb.FFT_chg.FFT3D(phi, -1);

    // normalization
    for(int j=0; j<pwb.ngmc; j++)
    {
        gradphi_x[j] /= pwb.ngmc;
        gradphi_y[j] /= pwb.ngmc;
        gradphi_z[j] /= pwb.ngmc;

        // gradphi_x[j].imag() /= pwb.ngmc;
        // gradphi_y[j].imag() /= pwb.ngmc;
        // gradphi_z[j].imag() /= pwb.ngmc;

        phi_work[j] /= pwb.ngmc;
        // phi_work[j].imag() /= pwb.ngmc;

        phi[j] /= pwb.ngmc;
        // phi[j].imag() /= pwb.ngmc;
    }

    // div(epsilon*grad phi) in kspace
    // add the kappa^2 contrib, -phi_work
    for(int ig=0; ig<pwb.ngmc; ig++)
    {
        // lp[ig].real() = ModuleBase::TWO_PI * 
        //     (gradphi_x[ig].real()*pwb.gcar[ig].x + gradphi_y[ig].real()*pwb.gcar[ig].y + gradphi_z[ig].real()*pwb.gcar[ig].z
        //     - phi_work[ig].real());
        
        // lp[ig].imag() = ModuleBase::TWO_PI * 
        //     (gradphi_x[ig].imag()*pwb.gcar[ig].x + gradphi_y[ig].imag()*pwb.gcar[ig].y + gradphi_z[ig].imag()*pwb.gcar[ig].z
        //     - phi_work[ig].imag());
        
        lp[ig] = ModuleBase::TWO_PI * 
            (gradphi_x[ig]*pwb.gcar[ig].x + gradphi_y[ig]*pwb.gcar[ig].y + gradphi_z[ig]*pwb.gcar[ig].z
            - phi_work[ig]);
    }
}

// Routine for solving poissons eqn using conjugate gradient (CG Method) ..
void H_Hartree_pw::minimize(
    const UnitCell &ucell,
    PW_Basis &pwb,
    double *d_eps,
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
    double eps_bar = 0;
    // r * r
    double r2 = 0;
    // precond loop parameter
    int i = 0;
    // malloc vectors
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

    for(int ig=0; ig<pwb.ngmc; ig++)
    {
        gg = pwb.get_NormG_cartesian(ig);
        gsqu[ig].real(1.0 / (gg * ucell.tpiba2)); // without kappa_2
        gsqu[ig].imag(0);
    }

    // init guess for phi
    for(int ig = 0; ig < pwb.ngmc; ig++)
    {
        // sol_phi = ?
        phi[ig].real(tot_N[ig].real() * gsqu[ig].real());
        phi[ig].imag(tot_N[ig].imag() * gsqu[ig].real()); 
    }
    
    // call Leps to calculate div ( epsilon * grad ) phi - kappa^2 * phi
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
    // calculate Lp: pwb.ngmc, complex-double

    // the residue
    // r = A*phi + (chtot + N)
    // pwb.ngmc = ?
    for(int ig = 0; ig < pwb.ngmc; ig++)
    {
        resid[ig].real(lp[ig].real() + tot_N[ig].real());
        resid[ig].imag(lp[ig].imag() + tot_N[ig].imag());
    }
    // precondition of the residue, z = invLr
    for(int ig = 0; ig < pwb.ngmc; ig++)
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
    for(int ig = 0; ig < pwb.ngmc; ig++)
    {
        d[ig].real(z[ig].real());
        d[ig].imag(z[ig].imag());
    }

    // CG Loop
    while(count < 20000 && sqrt(r2) > ediffsol && sqrt(rinvLr) > 1e-10)
    {
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
        // calculate alpha
        alpha = -rinvLr / Diago_CG::ddot_real(pwb.ngmc, d, lp);
        // update phi
        for(int ig = 0; ig < pwb.ngmc; ig++)
        {
            // phi[ig].real() += alpha * d[ig].real();
            // phi[ig].imag() += alpha * d[ig].imag();
            phi[ig] += alpha * d[ig];
        }
        // update resid
        for(int ig = 0; ig < pwb.ngmc; ig++)
        {
            // resid[ig].real() += alpha * lp[ig].real();
            // resid[ig].imag() += alpha * lp[ig].imag();
            resid[ig] += alpha * lp[ig];
        }
        // precond one more time..
        for(int ig = 0; ig < pwb.ngmc; ig++)
        {
            // z[ig].real() = gsqu[ig].real() * resid[ig].real();
            // z[ig].imag() = gsqu[ig].real() * resid[ig].imag();
            z[ig] = gsqu[ig] * resid[ig];
        }
        // calculate beta
        beta = 1.0 / rinvLr;
        rinvLr = Diago_CG::ddot_real(pwb.ngmc, resid, z);
        beta *= rinvLr;
        // update d
        for(int ig = 0; ig < pwb.ngmc; ig++)
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
    // output: num of cg loop
    ncgsol = count;

    // save phi for the intial guess for the next iteration
    // sol_phi
    // for(int ig = 0; ig < pwb.ngmc; ig++)
    // {
    //     sol_phi[ig].real(phi[ig].real());
    //     sol_phi[ig].imag(phi[ig].imag());
    // }
    // multiply by e/epsilon_0 and divide by volume for CHTOT
    
    // TODO:
    for(int ig = 0; ig < pwb.ngmc; ig++)
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

