#include "./kedf_vw.h"
#include <iostream>
#include "../src_parallel/parallel_reduce.h"

void KEDF_vW::set_para(int nx, double dV)
{
    this->nx = nx;
    this->dV = dV;
    if (this->potential != NULL)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            delete[] this->potential[is];
        }
        delete[] this->potential;
    } 
    this->potential = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->potential[is] = new double[this->nx];
        ModuleBase::GlobalFunc::ZEROS(this->potential[is], this->nx);
    }
}

// 
// EvW = -1/2 \int{sqrt(rho) nabla^2 sqrt(rho)}
// 
double KEDF_vW::get_energy(double **pphi, ModulePW::PW_Basis *pw_rho)
{
    double **LapPhi = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) LapPhi[is] = new double[pw_rho->nrxx];
    this->laplacianPhi(pphi, LapPhi, pw_rho);

    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pphi[0][ir] * LapPhi[0][ir];
        }
        energy *= this->dV * 0.5;
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx; ++ir)
            {
                energy += 2 * pphi[is][ir] * LapPhi[is][ir];
            }
        }
        energy *= 0.5 * this->dV * 0.5;
    }
    this->vWenergy = energy;
    Parallel_Reduce::reduce_double_all(this->vWenergy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] LapPhi[is];
    }
    delete[] LapPhi;

    return energy;
}

double KEDF_vW::get_energy_density(double **pphi, int is, int ir, ModulePW::PW_Basis *pw_rho)
{
    double **LapPhi = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) LapPhi[is] = new double[pw_rho->nrxx];
    this->laplacianPhi(pphi, LapPhi, pw_rho);

    double energyDen = 0.; // in Ry
    energyDen = 0.5 * pphi[is][ir] * LapPhi[is][ir];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] LapPhi[is];
    }
    delete[] LapPhi;
    return energyDen;
}

// 
// Vtf = delta Etf/delta rho = 5/3 * cTF * rho^{2/3}
// 
void KEDF_vW::vW_potential(const double * const *pphi, ModulePW::PW_Basis *pw_rho)
{
    ModuleBase::timer::tick("KEDF_vW", "vw_potential");

    double **LapPhi = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) LapPhi[is] = new double[pw_rho->nrxx];
    this->laplacianPhi(pphi, LapPhi, pw_rho);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            this->potential[is][ir] = 0.5 * LapPhi[is][ir] / pphi[is][ir];
        }
    }

    // calculate energy
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pphi[0][ir] * LapPhi[0][ir];
        }
        energy *= this->dV * 0.5;
    }
    else if (GlobalV::NSPIN == 2)
    {
        for (int is = 0; is < GlobalV::NSPIN; ++is)
        {
            for (int ir = 0; ir < this->nx; ++ir)
            {
                energy += 2 * pphi[is][ir] * LapPhi[is][ir];
            }
        }
        energy *= 0.5 * this->dV * 0.5;
    }
    this->vWenergy = energy;
    Parallel_Reduce::reduce_double_all(this->vWenergy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] LapPhi[is];
    }
    delete[] LapPhi;

    ModuleBase::timer::tick("KEDF_vW", "vw_potential");
}

void KEDF_vW::get_stress(double cellVol, double inpt_vWenergy)
{
    // waiting ...
}

void KEDF_vW::laplacianPhi(const double * const * pphi, double **rLapPhi, ModulePW::PW_Basis *pw_rho)
{
    std::complex<double> **recipPhi = new std::complex<double> *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipPhi[is] = new std::complex<double>[pw_rho->npw];

        pw_rho->real2recip(pphi[is], recipPhi[is]);
        for (int ik = 0; ik < pw_rho->npw; ++ik)
        {
            recipPhi[is][ik] *= pw_rho->gg[ik] * pw_rho->tpiba2 * 2; // mutiply by 2 to convert Hartree to Ry
        }
        pw_rho->recip2real(recipPhi[is], rLapPhi[is]);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] recipPhi[is];
    }
    delete[] recipPhi;
}