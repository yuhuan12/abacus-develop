#include "esolver_fp.h"
#include "./opt/opt_TN.hpp"
#include "./opt/opt_DCsrch.h"
#include "../module_psi/psi.h"

#include "./kedf_tf.h"
namespace ModuleESolver
{
class ESolver_OF: public ESolver_FP
{
// =========== TO DO LIST =============
// PARALLEL SUPPORT
// MORE KEDF
// MD TEST
// SPIN POLARISE
// REPLACE LINE SEARCH
// ADD UNIT TEST
public:

    // psi::Psi<double> *ppsi; 
    psi::Psi<double> ppsi; 
    // ElecState *p_es; 

    ESolver_OF()
    {
        // this->p_es = new ElecState_PW();
        // this->phamilt = new Hamilt_PW();
        this->classname = "ESolver_OF";

        // this->pdirect = new *double[1];
        // this->pphi = new *double[1];
        // this->pdLdphi = new *double[1];

        // this->theta = new double[1];
        // this->mu = new double[1];
        // this->normdLdphi = new double[1];
        this->task = new char[60];
    }

    ~ESolver_OF()
    {
        // delete this->p_es;
        // delete this->ppsi;
        if (this->pdirect != NULL)
        {
            for (int i = 0; i < GlobalV::NSPIN; ++i)
            {
                delete[] this->pdirect[i];
            }
            delete[] this->pdirect;
        } 
        if (this->pphi != NULL) delete[] this->pphi;
        if (this->pdLdphi != NULL)
        {
            for (int i = 0; i < GlobalV::NSPIN; ++i)
            {
                delete[] this->pdLdphi[i];
            }
            delete[] this->pdLdphi;
        } 
        if (this->pdEdphi != NULL)
        {
            for (int i = 0; i < GlobalV::NSPIN; ++i)
            {
                delete[] this->pdEdphi[i];
            }
            delete[] this->pdEdphi;
        } 

        if (this->nelec != NULL) delete[] this->nelec;
        if (this->theta != NULL) delete[] this->theta;
        if (this->mu != NULL) delete[] this->mu;
        if (this->task != NULL) delete[] this->task;
        if (this->opt_cg_mag != NULL) delete[] this->opt_cg_mag;
    }

    // Basis_PW basis_pw;
    virtual void Init(Input &inp, UnitCell_pseudo &cell) override;
    virtual void Run(int istep, UnitCell_pseudo& cell) override;

    virtual void cal_Energy(energy &en) override;
    virtual void cal_Force(ModuleBase::matrix &force) override;
    virtual void cal_Stress(ModuleBase::matrix &stress) override;

    virtual int getniter() override {
        return this->iter;
    }

private:
    KEDF_TF tf;

    Opt_CG opt_cg;
    Opt_TN opt_tn;
    Opt_DCsrch opt_dcsrch;

    Opt_CG *opt_cg_mag;

    // form Input
    string of_kinetic = "wt";   // Kinetic energy functional, such as TF, VW, WT
    string of_method = "tn";    // optimization method, include cg1, cg2, tn (default), bfgs
    string of_conv = "energy";  // select the convergence criterion, potential, energy (default), or both
    double of_tole = 2e-6;      // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
    double of_tolp = 1e-5;      // tolerance of potential for determining the convergence, default=1e-5 in a.u.
    double *nelec = NULL;              // number of electrons
    // double *nelecspin = NULL;   // number of spin up and spin down electrons
    int maxIter = 50;           // scf_nmax

    int iter = 0;
    int nrxx = 0; // PWBASIS
    double dV = 0; // CELL
    int maxDCsrch = 200; // max no. of line search

    // used in density optimization
    double **pdirect = NULL;
    double *theta = NULL;
    // double theta = 0.2;
    double **pdEdphi = NULL; // dE/dphi
    double **pdLdphi = NULL; // dL/dphi
    double **pphi = NULL; // pphi[i] = ppsi.get_pointer(i), which will be freed in ~Psi().
    char *task = NULL; // used in line search
    double *mu = NULL; // chemical potential
    double normdLdphi = 100.;
    int tnSpinFlag = -1; // spin flag used in calV, which will be called by opt_tn

    // used in conv check
    bool conv = false;
    double energy_llast = 0;
    double energy_last = 0;
    double energy_current = 0;

    int flag = -1; // flag of TN

    void updateV(); // waiting
    void solveV();
    void getNextDirect();
    void updateRho();
    bool checkExit();

    void calV(double *ptempPhi, double *rdLdphi);
    void caldEdtheta(double **ptempPhi, double **ptempRho, double *ptheta, double *rdEdtheta);
    double cal_mu(double *pphi, double *pdEdphi, double nelec);
    void printInfo();
    double inner_product(double *pa, double *pb, int length, double dV=1)
    {
        double innerproduct = 0.;
        for (int i = 0; i < length; ++i) innerproduct += pa[i] * pb[i];
        innerproduct *= dV;
        return innerproduct;
    }
};
}
