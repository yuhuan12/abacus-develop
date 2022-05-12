#include "esolver_fp.h"
#include "./opt/opt_TN.hpp"
#include "./opt/opt_DCsrch.h"
#include "../module_psi/psi.h"

#include "./kedf_tf.h"
namespace ModuleESolver
{
class ESolver_OF: public ESolver_FP
{
public:

    // psi::Psi<double> *ppsi; 
    psi::Psi<double> ppsi; 
    // ElecState *p_es; 

    ESolver_OF()
    {
        // this->p_es = new ElecState_PW();
        // this->phamilt = new Hamilt_PW();
        this->classname = "ESolver_OF";
        this->pdirect = new double[1];
        this->pphi = new double[1];
        this->pdLdphi = new double[1];
        this->task = new char[60];
    }

    ~ESolver_OF()
    {
        // delete this->p_es;
        // delete this->ppsi;
        if (this->pdirect != NULL) delete[] this->pdirect;
        if (this->pdLdphi != NULL) delete[] this->pdLdphi;
        if (this->task != NULL) delete[] this->task;
    }

    // Basis_PW basis_pw;
    virtual void Init(Input &inp, UnitCell_pseudo &cell) override;
    virtual void Run(int istep, UnitCell_pseudo& cell) override;

    virtual void cal_Energy(energy &en) override;  // waiting
    virtual void cal_Force(ModuleBase::matrix &force) override {}; // waiting
    virtual void cal_Stress(ModuleBase::matrix &stress) override {}; // waiting

    virtual int getniter() override {
        return this->iter;
    }

private:
    KEDF_TF tf;

    Opt_CG opt_cg;
    Opt_TN opt_tn;
    Opt_DCsrch opt_dcsrch;

    // form Input
    string of_kinetic = "wt";   // Kinetic energy functional, such as TF, VW, WT
    string of_method = "tn";    // optimization method, include cg1, cg2, tn (default), bfgs
    string of_conv = "energy";  // select the convergence criterion, potential, energy (default), or both
    double of_tole = 2e-6;      // tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry
    double of_tolp = 1e-5;      // tolerance of potential for determining the convergence, default=1e-5 in a.u.
    int nelec = 0;              // number of electrons
    int maxIter = 50;           // scf_nmax

    int iter = 0;
    int nrxx = 0; // PWBASIS
    double dV = 0; // CELL

    // used in density optimization
    double *pdirect = NULL;
    double theta = 0;
    double *pdLdphi = NULL; // dL/dphi
    double *pphi = NULL; // pphi = ppsi.get_pointer(0), it will be freed in ~Psi().
    char *task = NULL; // used in line search
    double mu = 0; // chemical potential

    // used in conv check
    bool conv = false;
    double energy_llast = 0;
    double energy_last = 0;
    double energy_current = 0;

    int flag = -1;

    void updateV(); // waiting
    void solveV();
    void updateRho();
    bool checkExit();

    void calV(double *ptempPhi, double *rdLdphi);
    double caldEdtheta(double *ptempPhi, double **ptempRho);
    double cal_mu(double *pphi, double *pdEdphi);
    double inner_product(double *pa, double *pb, int length, double dV=1)
    {
        double innerproduct = 0.;
        for (int i = 0; i < length; ++i) innerproduct += pa[i] * pb[i];
        innerproduct *= dV;
        return innerproduct;
    }
};
}