#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
#include "module_hamilt_pw/hamilt_pwdft/operator_pw/velocity_pw.h"
// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"


namespace ModuleESolver
{

    template<typename FPTYPE, typename Device = psi::DEVICE_CPU>
    class ESolver_KS_PW : public ESolver_KS<FPTYPE, Device>
    {
    public:
        ESolver_KS_PW();
        ~ESolver_KS_PW();
        void Init(Input& inp, UnitCell& cell) override;
        void cal_Energy(double& etot) override;
        void cal_Force(ModuleBase::matrix& force) override;
        void cal_Stress(ModuleBase::matrix& stress) override;
        virtual void hamilt2density(const int istep, const int iter, const FPTYPE ethr) override;
        virtual void hamilt2estates(const FPTYPE ethr) override;
        virtual void nscf() override;
        void postprocess() override;
        //calculate conductivities with Kubo-Greenwood formula
        void KG(const int nche_KG, const double fwhmin, const double wcut,
             const double dw_in, const double dt_in, ModuleBase::matrix& wg);
        void jjcorr_ks(const int ik, const int nt, const double dt, ModuleBase::matrix& wg, hamilt::Velocity& velop, 
                       double* ct11, double* ct12, double* ct22);

    protected:
        virtual void beforescf(const int istep) override;
        virtual void eachiterinit(const int istep, const int iter) override;
        virtual void updatepot(const int istep, const int iter) override;
        virtual void eachiterfinish(const int iter) override;
        virtual void afterscf(const int istep) override;
        virtual void othercalculation(const int istep)override;

        //temporary, this will be removed in the future;
        //Init Global class
        void Init_GlobalC(Input& inp, UnitCell& cell);
        //calculate conductivities from j-j correlation function
        void calcondw(const int nt,const double dt, const double fwhmin, const double wcut, const double dw_in, double *ct11, double *ct12, double *ct22);


    private:
        
        Device * ctx = {};
        psi::AbacusDevice_t device = {};
        psi::Psi<std::complex<FPTYPE>, Device>* kspw_psi = nullptr;
        psi::Psi<std::complex<double>, Device>* __kspw_psi = nullptr;
        using castmem_2d_d2h_op = psi::memory::cast_memory_op<std::complex<double>, std::complex<FPTYPE>, psi::DEVICE_CPU, Device>;
    };
}  // namespace ModuleESolver
#endif
