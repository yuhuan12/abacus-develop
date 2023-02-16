#ifndef ELECSTATELCAOTDDFT_H
#define ELECSTATELCAOTDDFT_H

#include "elecstate.h"
#include "elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"

namespace elecstate
{

class ElecStateLCAO_TDDFT : public ElecStateLCAO
{
  public:
    ElecStateLCAO_TDDFT(Charge* chg_in = nullptr,
                        const K_Vectors* klist_in = nullptr,
                        int nks_in = 1,
                        Local_Orbital_Charge* loc_in = nullptr,
                        LCAO_Hamilt* uhm_in = nullptr,
                        Local_Orbital_wfc* lowf_in = nullptr)
    {
        init_ks(chg_in, klist_in, nks_in);
        this->loc = loc_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
        this->classname = "ElecStateLCAO_TDDFT";
    }
    void psiToRho_td(const psi::Psi<std::complex<double>>& psi);
    void calculate_weights_td();
};

} // namespace elecstate

#endif