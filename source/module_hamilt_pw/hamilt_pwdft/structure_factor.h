#ifndef PLANEWAVE_H
#define PLANEWAVE_H

#include "module_cell/unitcell.h"
#include "module_base/complexmatrix.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_psi/psi.h"

using namespace std;

class Structure_Factor
{

public:
    Structure_Factor();
    ~Structure_Factor();
    void set(const int &nbspline_in);


//===============================================
// Part 4: G vectors in reciprocal FFT box
//===============================================
public:
    int nbspline;

	// structure factor (ntype, ngmc)
    ModuleBase::ComplexMatrix strucFac;
    void setup_structure_factor(UnitCell* Ucell,ModulePW::PW_Basis* rho_basis); 		// Calculate structure factors
    void bspline_sf(const int,UnitCell* Ucell,ModulePW::PW_Basis* rho_basis); //calculate structure factors through Cardinal B-spline interpolation
    void bsplinecoef(std::complex<double> *b1, std::complex<double> *b2, std::complex<double> *b3, 
                    const int nx, const int ny, const int nz, const int norder);


public:
	// phase of e^{-iG*tau_s}
    ModuleBase::ComplexMatrix eigts1; // dimension: [Ucell->nat, 2*this->ncx + 1] 
    ModuleBase::ComplexMatrix eigts2; // dimension: [Ucell->nat, 2*this->ncy + 1] 
    ModuleBase::ComplexMatrix eigts3; // dimension: [Ucell->nat, 2*this->ncz + 1]

    template <typename FPTYPE> std::complex<FPTYPE> * get_eigts1_data() const;
    template <typename FPTYPE> std::complex<FPTYPE> * get_eigts2_data() const;
    template <typename FPTYPE> std::complex<FPTYPE> * get_eigts3_data() const;

private:

    std::complex<float> * c_eigts1 = nullptr, * c_eigts2 = nullptr, * c_eigts3 = nullptr;
    std::complex<double> * z_eigts1 = nullptr, * z_eigts2 = nullptr, * z_eigts3 = nullptr;

};
#endif //PlaneWave class
