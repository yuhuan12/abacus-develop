#include "./kedf_tf.h"
#include <iostream>

// 
// Etf = cTF * \int{dr rho^{5/3}}
// 
double KEDF_TF::get_energy(double **prho)
{
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pow(prho[0][ir], 5./3.);
        }
        energy *= this->dV * this->cTF;
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx; ++ir)
            {
                energy += pow(2. * prho[is][ir], 5./3.);
            }
        }
        energy *= 0.5 * this->dV * this->cTF;
    }
    this->TFenergy = energy;
    return energy;
}

double KEDF_TF::get_energy_density(double **prho, int is, int ir)
{
    double energyDen = 0.; // in Ry
    energyDen = this->cTF * pow(prho[is][ir], 5./3.);
    return energyDen;
}

// 
// Vtf = delta Etf/delta rho = 5/3 * cTF * rho^{2/3}
// 
void KEDF_TF::get_potential(const double * const *prho)
{
    if (GlobalV::NSPIN  == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            this->potential[0][ir] = 5.0/3.0 * this->cTF * pow(prho[0][ir], 2./3.);
        }
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx; ++ir)
            {
                this->potential[is][ir] = 5.0/3.0 * this->cTF * pow(2. * prho[is][ir], 2./3.);
            }
        }
    }
}

void KEDF_TF::get_stress(double cellVol, double inpt_TFenergy)
{
    double temp = 0.;
    if (inpt_TFenergy == -1)
    {
        temp = 2. * this->TFenergy / (3. * cellVol);
    }
    else
    {
        temp = 2. * inpt_TFenergy / (3. * cellVol);
    }
// #ifdef __MPI
//     //
// #endif
    for (int i = 0; i < 3; ++i)
    {
        this->stress(i, i) = temp;
    }
}