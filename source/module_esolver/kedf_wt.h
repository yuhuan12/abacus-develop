#include <stdio.h>
#include <math.h>
#include "../module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_pw/pw_basis.h"

class KEDF_WT
{
public:
    KEDF_WT()
    {
        this->potential = NULL;
        this->kernel = NULL;
        this->stress.create(3,3);
    }
    ~KEDF_WT()
    {
        if (this->potential != NULL)
        {
            for (int is = 0; is < GlobalV::NSPIN; ++is)
            {
                delete[] this->potential[is];
            }
            delete[] this->potential;
        } 
        if (this->kernel != NULL) delete[] this->kernel;
    }

    void set_para(int nx, double dV, double nelec, double weightWT, double weightTF, double weightVW, ModulePW::PW_Basis *pw_rho);

    double get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho);
    double get_energy_density(double **prho, int is, int ir, ModulePW::PW_Basis *pw_rho);
    void WT_potential(const double * const * prho, ModulePW::PW_Basis *pw_rho);
    void get_stress(double cellVol, double inpt_vWenergy=-1);
    double WTkernel(double eta, double weightTF, double weightVW);
    void multiKernel(const double * const * prho, double **rkernelRho, double exponent, ModulePW::PW_Basis *pw_rho);
    // double getFermiWaveVector(double **prho);

    // int NSPIN = 1;
    int nx = 0;
    double dV = 0.;
    double kF = 0.;
    double alpha = (5. + sqrt(5.))/6.;
    double beta = (5. - sqrt(5.))/6.;
    double weightWT = 1.;
    const double cTF = 3.0/10.0 * pow(3*pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    double WTenergy = 0.;
    double **potential;
    double *kernel;
    ModuleBase::matrix stress;
};