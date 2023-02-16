#ifndef WRITE_PROJ_BAND_LCAO_H
#define WRITE_PROJ_BAND_LCAO_H
#include "module_psi/psi.h"
#include "module_elecstate/elecstate.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_cell/klist.h"
#include "module_cell/unitcell.h"
#include "module_orbital/ORB_read.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"

namespace ModuleIO
{
	void write_proj_band_lcao(const psi::Psi<double> *psid,
		const psi::Psi<std::complex<double>> *psi,
		LCAO_Hamilt& uhm,
		const elecstate::ElecState* pelec,
		const K_Vectors* kv,
		const UnitCell &ucell, 
		const LCAO_Orbitals &ORB,
		Grid_Driver &GridD);
}

#endif
