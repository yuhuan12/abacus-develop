#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_reduce.h"
#include "gint_k.h"
#include "module_basis/module_ao/ORB_read.h"
#include "grid_technique.h"
#include "module_base/ylm.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h" // mohan add 2021-01-30
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"

void Gint_k::allocate_pvdpR(void)
{
    ModuleBase::TITLE("Gint_k","allocate_pvpR");

    //xiaohui modify 2015-05-30
    // the number of matrix element <phi_0 | V | dphi_R> is GlobalC::GridT.nnrg.
    this->pvdpRx_reduced = new double*[GlobalV::NSPIN];
    this->pvdpRy_reduced = new double*[GlobalV::NSPIN];
    this->pvdpRz_reduced = new double*[GlobalV::NSPIN];
    for(int is =0;is<GlobalV::NSPIN;is++)
    {
        this->pvdpRx_reduced[is] = new double[GlobalC::GridT.nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvdpRx_reduced[is], GlobalC::GridT.nnrg);
        this->pvdpRy_reduced[is] = new double[GlobalC::GridT.nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvdpRy_reduced[is], GlobalC::GridT.nnrg);
        this->pvdpRz_reduced[is] = new double[GlobalC::GridT.nnrg];	
        ModuleBase::GlobalFunc::ZEROS( pvdpRz_reduced[is], GlobalC::GridT.nnrg);
    }

    ModuleBase::Memory::record("pvdpR_reduced", 3 * sizeof(double) * GlobalC::GridT.nnrg * GlobalV::NSPIN);

    return;
}

void Gint_k::destroy_pvdpR(void)
{
    ModuleBase::TITLE("Gint_k","destroy_pvpR");
    
    for(int is =0;is<GlobalV::NSPIN;is++) delete[] pvdpRx_reduced[is];
    delete[] pvdpRx_reduced;
    for(int is =0;is<GlobalV::NSPIN;is++) delete[] pvdpRy_reduced[is];
    delete[] pvdpRy_reduced;
    for(int is =0;is<GlobalV::NSPIN;is++) delete[] pvdpRz_reduced[is];
    delete[] pvdpRz_reduced;

    return;
}