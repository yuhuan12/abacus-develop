#include "./kedf_wt.h"
#include <iostream>
#include "../src_parallel/parallel_reduce.h"

void KEDF_WT::set_para(int nx, double dV, double nelec, double weightWT, double weightTF, double weightVW, ModulePW::PW_Basis *pw_rho)
{
    this->nx = nx;
    this->dV = dV;
    this->weightWT = weightWT;
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

    double rho0 = 1./(nx * dV) * nelec;
    this->kF = pow(3. * pow(ModuleBase::PI, 2) * rho0, 1./3.);

    if (this->kernel != NULL) delete[] this->kernel;
    this->kernel = new double[pw_rho->npw];
    double eta = 0.;
    for (int ip = 0; ip < pw_rho->npw; ++ip)
    {
        eta = sqrt(pw_rho->gg[ip]) * pw_rho->tpiba / this->kF / 2.;
        this->kernel[ip] = this->WTkernel(eta, weightTF, weightVW);
    }
}

// 
// Ewt = Ctf * \int{rho^alpha * W(r - r') * rho^beta drdr'} + T_vW + T_TF, T_vW and T_TF will be added in ESolver_OF
// 
double KEDF_WT::get_energy(const double * const * prho, ModulePW::PW_Basis *pw_rho)
{
    double **kernelRhoBeta = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoBeta, this->beta, pw_rho);
    // std::complex<double> **recipkernelRhoBeta = new std::complex<double> *[GlobalV::NSPIN];
    // for (int is = 0; is < GlobalV::NSPIN; ++is)
    // {
    //     kernelRhoBeta[is] = new double[pw_rho->nrxx];
    //     recipkernelRhoBeta[is] = new std::complex<double>[pw_rho->npw];
    //     for (int ir = 0; ir < this->nx; ++ir)
    //     {
    //         kernelRhoBeta[is][ir] = pow(prho[is][ir], this->beta);
    //     }
    //     pw_rho->real2recip(kernelRhoBeta, recipkernelRhoBeta);
    //     for (int ip = 0; ip < pw_rho->npw; ++ip)
    //     {
    //         recipkernelRhoBeta[is][ip] *= this->kernel[ip];
    //     }
    //     pw_rho->recip2real(recipkernelRhoBeta, kernelRhoBeta);
    // } 

    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pow(prho[0][ir], this->alpha) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV * this->cTF * this->weightWT;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // for (int is = 0; is < GlobalV::NSPIN; ++is)
        // {
        //     for (int ir = 0; ir < this->nx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV * 0.5;
    }
    this->WTenergy = energy;
    Parallel_Reduce::reduce_double_all(this->WTenergy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;

    return energy;
}

double KEDF_WT::get_energy_density(double **prho, int is, int ir, ModulePW::PW_Basis *pw_rho)
{
    double **kernelRhoBeta = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoBeta, this->beta, pw_rho);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
    }
    delete[] kernelRhoBeta;
    return this->cTF * this->weightWT * (prho[is][ir], this->alpha) * kernelRhoBeta[is][ir];
}


// 
// Vtf = delta Etf/delta rho = 5/3 * cTF * rho^{2/3}
// 
void KEDF_WT::WT_potential(const double * const *prho, ModulePW::PW_Basis *pw_rho)
{
    ModuleBase::timer::tick("KEDF_WT", "wt_potential");

    double **kernelRhoBeta = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoBeta[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoBeta, this->beta, pw_rho);

    double **kernelRhoAlpha = new double* [GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is) kernelRhoAlpha[is] = new double[pw_rho->nrxx];
    this->multiKernel(prho, kernelRhoAlpha, this->alpha, pw_rho);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            this->potential[is][ir] = this->cTF * this->weightWT * 
                                    (this->alpha * pow(prho[is][ir], this->alpha-1.) * kernelRhoBeta[is][ir]
                                    + this->beta *pow(prho[is][ir], this->beta-1.) * kernelRhoAlpha[is][ir]);
        }
    }

    // calculate energy
    double energy = 0.; // in Ry
    if (GlobalV::NSPIN == 1)
    {
        for (int ir = 0; ir < this->nx; ++ir)
        {
            energy += pow(prho[0][ir], this->alpha) * kernelRhoBeta[0][ir];
        }
        energy *= this->dV * this->cTF * this->weightWT;
    }
    else if (GlobalV::NSPIN == 2)
    {
        // for (int is = 0; is < GlobalV::NSPIN; ++is)
        // {
        //     for (int ir = 0; ir < this->nx; ++ir)
        //     {
        //         energy += 2 * pphi[is][ir] * LapPhi[is][ir];
        //     }
        // }
        // energy *= 0.5 * this->dV * 0.5;
    }
    this->WTenergy = energy;
    Parallel_Reduce::reduce_double_all(this->WTenergy);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] kernelRhoBeta[is];
        delete[] kernelRhoAlpha[is];
    }
    delete[] kernelRhoBeta;
    delete[] kernelRhoAlpha;
    ModuleBase::timer::tick("KEDF_WT", "wt_potential");
}

void KEDF_WT::get_stress(double cellVol, double inpt_vWenergy)
{
    // waiting ...
}

double KEDF_WT::WTkernel(double eta, double weightTF, double weightVW)
{
    if (eta < 0.) 
    {
        return 0.;
    }
    // limit for small eta
    else if (eta < 1e-10)
    {
        return 1. - weightTF + eta * eta * (1./3. - 3. * weightVW);
    }
    // around the singularity
    else if (abs(eta - 1.) < 1e-10)
    {
        return 2. - weightTF - 3. * weightVW + 20. * (eta - 1);
    }
    // Taylor expansion for high eta
    else if (eta > 3.65)
    {
        double eta2 = eta * eta;
        double invEta2 = 1. / eta2;
        double LindG = 3. * (1. - weightVW) * eta2 
                        -weightTF-0.6 
                        + invEta2 * (-0.13714285714285712 
                        + invEta2 * (-6.39999999999999875E-2
                        + invEta2 * (-3.77825602968460128E-2
                        + invEta2 * (-2.51824061652633074E-2
                        + invEta2 * (-1.80879839616166146E-2
                        + invEta2 * (-1.36715733124818332E-2
                        + invEta2 * (-1.07236045520990083E-2
                        + invEta2 * (-8.65192783339199453E-3 
                        + invEta2 * (-7.1372762502456763E-3
                        + invEta2 * (-5.9945117538835746E-3
                        + invEta2 * (-5.10997527675418131E-3 
                        + invEta2 * (-4.41060829979912465E-3 
                        + invEta2 * (-3.84763737842981233E-3 
                        + invEta2 * (-3.38745061493813488E-3 
                        + invEta2 * (-3.00624946457977689E-3)))))))))))))));
        return LindG;
    }
    else
    {
        return 1. / (0.5 + 0.25 * (1. - eta * eta) / eta * log((1 + eta)/abs(1 - eta)))
                 - 3. * weightVW * eta * eta - weightTF;
    }
}


void KEDF_WT::multiKernel(const double * const * prho, double **rkernelRho, double exponent, ModulePW::PW_Basis *pw_rho)
{
    std::complex<double> **recipkernelRho = new std::complex<double> *[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        recipkernelRho[is] = new std::complex<double>[pw_rho->npw];
        for (int ir = 0; ir < this->nx; ++ir)
        {
            rkernelRho[is][ir] = pow(prho[is][ir], exponent);
        }
        pw_rho->real2recip(rkernelRho[is], recipkernelRho[is]);
        for (int ip = 0; ip < pw_rho->npw; ++ip)
        {
            recipkernelRho[is][ip] *= this->kernel[ip];
        }
        pw_rho->recip2real(recipkernelRho[is], rkernelRho[is]);
    }

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        delete[] recipkernelRho[is];
    }
    delete[] recipkernelRho;
}