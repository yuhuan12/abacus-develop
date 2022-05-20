#include <stdio.h>
#include <math.h>
#include "../module_base/global_function.h"

class KEDF_TF
{
public:
    KEDF_TF()
    {
        this->potential = NULL;
    }
    ~KEDF_TF()
    {
        if (this->potential != NULL) delete[] this->potential;
    }

    void set_para(int nx, double dV)
    {
        this->nx = nx;
        this->dV = dV;
        delete[] this->potential;
        this->potential = new double[this->nx];
        ModuleBase::GlobalFunc::ZEROS(this->potential, this->nx);
    }
    double get_energy(double **prho);
    void get_potential(const double * const * prho);

    int NSPIN = 1;
    int nx = 0;
    double dV = 0.;
    // const double cTF = 0.; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    const double cTF = 3.0/10.0 * pow(3*pow(M_PI, 2.0), 2.0/3.0) * 2; // 10/3*(3*pi^2)^{2/3}, multiply by 2 to convert unit from Hartree to Ry, finally in Ry*Bohr^(-2)
    double *potential;
};