#include <stdio.h>
#include <math.h>
#include "../module_base/timer.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"

class KEDF_TF
{
public:
    KEDF_TF()
    {
        this->potential = NULL;
        this->stress.create(3,3);
    }
    ~KEDF_TF()
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

    double get_energy(double **prho);
    double get_energy_density(double **prho, int is, int ir);
    void tf_potential(const double * const * prho);
    void get_stress(double cellVol, double inpt_TFenergy=-1);


    // int NSPIN = 1;
    int nx = 0;
    double dV = 0.;
    // const double cTF = 0.; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double cTF = 3.0/10.0 * pow(3*pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    double TFenergy = 0.;
    double **potential;
    ModuleBase::matrix stress;
};