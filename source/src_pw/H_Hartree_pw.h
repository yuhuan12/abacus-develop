#ifndef H_HARTREE_PW_H
#define H_HARTREE_PW_H

#include "tools.h"
#include "../module_cell/unitcell.h"
#include "pw_basis.h"
#include "use_fft.h"

class H_Hartree_pw 
{
	public:

    H_Hartree_pw();
    ~H_Hartree_pw();

	// the Hartree energy
    static double hartree_energy;

	// compute the Hartree energy
    static ModuleBase::matrix v_hartree(
		const UnitCell &cell, 
		PW_Basis &pwb, 
		const int &nspin,
		const double*const*const rho);

	void gauss_charge(
		const UnitCell &cell, 
		PW_Basis &pwb, 
		complex<double> *N,
		const int flag);
	
	int get_Z(string str);

	void V_correction(
		const UnitCell &cell,
		PW_Basis &pwb, 
		double **rho);

	void cast_C2R(complex<double> *src, double* dst, int dim);

	void Leps(
		const UnitCell &ucell,
		PW_Basis &pwb,
		complex<double> *phi,
		double *epsilon, // epsilon from shapefunc
		complex<double> *gradphi_x,
		complex<double> *gradphi_y,
		complex<double> *gradphi_z,
		complex<double> *phi_work,
		complex<double> *lp // output
	);
	
	void minimize(
		const UnitCell &ucell,
		PW_Basis &pwb,
		double *d_eps,
		const complex<double>* tot_N,
		complex<double> *phi,
		int &ncgsol);

	private:


};

#endif //Hartree energy
