#include "exx_abfs.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h"
#include "module_base/mathzone.h"
#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_base/container_operator.h"

#include "src_ri/test_code/matrix-test.h"
#include "src_ri/test_code/complexmatrix-test.h"
#include "src_ri/test_code/exx_lcao-test.h"
#include "src_ri/test_code/exx_abfs-dm-test.h"

void Exx_Abfs::DM::cal_DM(
	const ModuleBase::matrix& wg,
	const std::set<std::pair<size_t,size_t>> &atom_pairs,
    const std::vector<Abfs::Vector3_Order<int>>& Born_von_Karman_boxes,
    std::complex<double>*** wfc_k_grid)
{
	ModuleBase::TITLE("Exx_Abfs::DM::cal_DM");
	
	cal_DMk_mixing( GlobalC::CHR_MIX, wg, atom_pairs, wfc_k_grid);

	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];

		for( const ModuleBase::Vector3<int> &box : Born_von_Karman_boxes )
		{
			DMr[iat1][iat2][box] = std::vector<ModuleBase::matrix>( GlobalV::NSPIN, {GlobalC::ucell.atoms[it1].nw,GlobalC::ucell.atoms[it2].nw} );
			for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
			{
				DMr[iat1][iat2][box][GlobalC::kv.isk[ik]] += ( DMk[iat1][iat2][ik] * exp( -ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT* (GlobalC::kv.kvec_c[ik]* (box*GlobalC::ucell.latvec)) ) ).real();
			}
		}
	}
}



void Exx_Abfs::DM::cal_DMk_mixing(
	const Charge_Mixing &charge,
	const ModuleBase::matrix& wg,
    const std::set<std::pair<size_t, size_t>>& atom_pairs,
     std::complex<double>*** wfc_k_grid)
{
	ModuleBase::TITLE("Exx_Abfs::DM::cal_DMk_mixing");

	if(flag_mix)
	{
		if ( charge.get_mixing_mode() == "plain" )
		{
			plain_mixing( charge, wg, atom_pairs,wfc_k_grid );
		}
		else if ( charge.get_mixing_mode() == "pulay" )
		{
			pulay_mixing( charge, wg, atom_pairs, wfc_k_grid );
		}
		else
		{
			throw std::invalid_argument("mixing density matrix error. In "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
		}
	}
	else
	{
		DMk = cal_DMk_raw(atom_pairs, wg, wfc_k_grid);
//		DMk = Exx_Abfs_DM_Test::cal_DMk_raw_readfile(atom_pairs);			// Peize Lin test 2018-03-27

		#if TEST_EXX_LCAO==1
			static int istep=0;
			ofs_matrixes("DMk_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+".dat(@Exx_Abfs::DM::cal_DMk_mixing)",DMk);
		#elif TEST_EXX_LCAO==-1
			#error
		#endif
	}
}



std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>> Exx_Abfs::DM::cal_DMk_raw(
	const std::set<std::pair<size_t,size_t>> &atom_pairs, 
	const ModuleBase::matrix& wg, 
	std::complex<double>*** wfc_k_grid ) const
{
	ModuleBase::TITLE("Exx_Abfs::DM::cal_DMk_raw");

	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	
	std::map<size_t,std::map<size_t,std::vector<ModuleBase::ComplexMatrix>>> DMk_raw;
	for( const std::pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = GlobalC::ucell.iat2it[iat1];
		const size_t it2 = GlobalC::ucell.iat2it[iat2];
		const size_t ia1 = GlobalC::ucell.iat2ia[iat1];
		const size_t ia2 = GlobalC::ucell.iat2ia[iat2];

		DMk_raw[iat1][iat2] = std::vector<ModuleBase::ComplexMatrix>( GlobalC::kv.nks, {GlobalC::ucell.atoms[it1].nw,GlobalC::ucell.atoms[it2].nw} );
		for( size_t ik=0; ik!=GlobalC::kv.nks; ++ik )
		{
			for( size_t iw1=0; iw1!=GlobalC::ucell.atoms[it1].nw; ++iw1 )
			{
				for( size_t iw2=0; iw2!=GlobalC::ucell.atoms[it2].nw; ++iw2 )
				{
					for( size_t ib=0; ib!=GlobalV::NBANDS; ++ib )
					{
						if( GlobalV::GAMMA_ONLY_LOCAL )
						{
							ModuleBase::WARNING_QUIT("Exx_Abfs::DM::cal_DMk_raw","need to update GlobalC::LOWF.WFC_GAMMA");
						}
						else
						{
							DMk_raw[iat1][iat2][ik](iw1,iw2) += wg(ik,ib) 
								* wfc_k_grid[ik][ib][GlobalC::ucell.itiaiw2iwt(it1,ia1,iw1)] 
								* conj(wfc_k_grid[ik][ib][GlobalC::ucell.itiaiw2iwt(it2,ia2,iw2)]);
						}
					}
				}
			}
			DMk_raw[iat1][iat2][ik] *= SPIN_multiple;
		}
	}
	
	#if TEST_EXX_LCAO==1
		static int istep=0;
		ofs_matrixes("DMk_raw_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+".dat(@Exx_Abfs::DM::cal_DMk_raw)",DMk_raw);
	#elif TEST_EXX_LCAO==-1
		#error
	#endif
	return DMk_raw;
}



void Exx_Abfs::DM::plain_mixing(
	const Charge_Mixing &charge,
	const ModuleBase::matrix& wg,
    const std::set<std::pair<size_t, size_t>>& atom_pairs,
    complex<double>*** wfc_k_grid)
{
	ModuleBase::TITLE("Exx_Abfs::DM::plain_mixing");

	if(DMk.empty())
		DMk = cal_DMk_raw(atom_pairs, wg, wfc_k_grid);
	else
		DMk = charge.get_mixing_beta() * cal_DMk_raw(atom_pairs, wg, wfc_k_grid) + (1-charge.get_mixing_beta()) * DMk;
}



void Exx_Abfs::DM::pulay_mixing(
	const Charge_Mixing &charge,
	const ModuleBase::matrix& wg,
    const std::set<std::pair<size_t, size_t>>& atom_pairs,
    complex<double>*** wfc_k_grid)
{
	if( 1==charge.get_totstep() )
	{
		DMk_pulay_seq.clear();
	}
	
	DMk_pulay_seq.push_back( charge.get_mixing_beta() * cal_DMk_raw(atom_pairs, wg, wfc_k_grid) + (1-charge.get_mixing_beta()) * DMk );
	if( charge.get_totstep() > charge.get_rstep() )
		DMk_pulay_seq.pop_front();
	
	if( 1==charge.get_totstep() )
	{
		DMk = DMk_pulay_seq.front();
	}
	else
	{
		const int alpha_size = DMk_pulay_seq.size()-1;
		const int alpha_end = charge.get_idstep();
		const int alpha_begin = alpha_end - alpha_size;
		auto alpha_index = [&](const int i){ return ((alpha_begin+i)%charge.get_dstep()+charge.get_dstep())%charge.get_dstep(); };
		
		DMk = (1+charge.get_alpha()[alpha_index(alpha_size-1)]) * DMk_pulay_seq.back();
		for( size_t i=1; i<DMk_pulay_seq.size()-1; ++i )
			DMk = DMk + ( charge.get_alpha()[alpha_index(i-1)] - charge.get_alpha()[alpha_index(i)] ) * DMk_pulay_seq[i];
		DMk = DMk - charge.get_alpha()[alpha_index(0)] * DMk_pulay_seq.front();
				
		#if TEST_EXX_LCAO==1
		{
			std::cout<<"charge.get_alpha()"<<std::endl;
			std::cout<<alpha_begin<<"\t"<<alpha_end<<"\t"<<alpha_size<<std::endl;
			std::cout<<charge.get_alpha()[alpha_index(alpha_size-1)]<<std::endl;
			for( size_t i=1; i<DMk_pulay_seq.size()-1; ++i )
				std::cout<<charge.get_alpha()[alpha_index(i-1)]<<"\t"<<charge.get_alpha()[alpha_index(i)]<<std::endl;
			std::cout<<charge.get_alpha()[alpha_index(0)]<<std::endl;
		}
		#elif TEST_EXX_LCAO==-1
			#error
		#endif
	}
	
	#if TEST_EXX_LCAO==1
	{
		std::cout<<"charge.alpha_all"<<std::endl;
		for( size_t i=0; i!=charge.get_dstep(); ++i )
			std::cout<<charge.get_alpha()[i]<<"\t";
		std::cout<<std::endl;
	}
	#elif TEST_EXX_LCAO==-1
		#error
	#endif	
	
	#if TEST_EXX_LCAO==1
//	{
//		static int istep=0;
//		for( size_t i=0; i!=DMk_pulay_seq.size(); ++i )
//			ofs_matrixes( "DMk_pulay_seq_"+ModuleBase::GlobalFunc::TO_STRING(istep)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), DMk_pulay_seq[i] );
//		++istep;
//	}
	#elif TEST_EXX_LCAO==-1
		#error
	#endif		
}
