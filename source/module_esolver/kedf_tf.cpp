#include "./kedf_tf.h"
#include <iostream>

// 
// Etf = cTF * \int{dr rho^{5/3}}
// 
double KEDF_TF::get_energy(double **prho)
{
    double energy = 0.; // in Ry
    if (NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pow(prho[0][ir], 5./3.);
        }
        energy *= this->dV * this->cTF;
    }
    else {}
    return energy;
}

// 
// Vtf = delta Etf/delta rho = 5/3 * cTF * rho^{2/3}
// 
void KEDF_TF::get_potential(const double * const *prho)
{
    if (NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            this->potential[ir] = 5.0/3.0 * this->cTF * pow(prho[0][ir], 2./3.);
        }
    }
}