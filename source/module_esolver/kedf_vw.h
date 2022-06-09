#include <stdio.h>
#include <math.h>
#include "../module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_pw/pw_basis.h"

class KEDF_vW
{
public:
    KEDF_vW()
    {
        this->potential = NULL;
        this->stress.create(3,3);
    }
    ~KEDF_vW()
    {
        if (this->potential != NULL)
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                delete[] this->potential[is];
            }
            delete[] this->potential;
        } 
    }

    void set_para(int nx, double dV);

    double get_energy(double **pphi, ModulePW::PW_Basis *pw_rho);
    double get_energy_density(double **pphi, int is, int ir, ModulePW::PW_Basis *pw_rho);
    void vW_potential(const double * const * pphi, ModulePW::PW_Basis *pw_rho);
    void get_stress(double cellVol, double inpt_vWenergy=-1);
    void laplacianPhi(const double * const * pphi, double **rLapPhi, ModulePW::PW_Basis *pw_rho);

    // int NSPIN = 1;
    int nx = 0;
    double dV = 0.;
    double vWenergy = 0.;
    double **potential;
    ModuleBase::matrix stress;
};