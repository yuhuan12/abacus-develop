#include "unitcell_pseudo.h"

#ifdef __LCAO
#include "../module_orbital/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif
//#include "../src_pw/global.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

UnitCell_pseudo::UnitCell_pseudo()
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("unitcell_pseudo","Constructor");
	set_atom_flag = false;
}

UnitCell_pseudo::~UnitCell_pseudo()
{
	if(set_atom_flag)
	{
		delete[] atoms;
	}
}

//==============================================================
//Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell_pseudo::setup_cell(
		const std::string &s_pseudopot_dir,
		output &outp,  
		const std::string &fn,
		std::ofstream &log)
{
	ModuleBase::TITLE("UnitCell_pseudo","setup_cell");	
	// (1) init mag
	assert(ntype>0);
#ifndef __CMD
	if(!input_mag)
	{
		delete[] magnet.start_magnetization;
	magnet.start_magnetization = new double[this->ntype];
	}	
#endif

	// (2) init *Atom class array.
	this->atoms = new Atom[this->ntype]; // atom species.
	this->set_atom_flag = true;


	bool ok = true;
	bool ok2 = true;

	// (3) read in atom information
	if(GlobalV::MY_RANK == 0)
	{
		// open "atom_unitcell" file.
		std::ifstream ifa(fn.c_str(), ios::in);
		if (!ifa)
		{
			GlobalV::ofs_warning << fn;
			ok = false;
		}

		if(ok)
		{

			GlobalV::ofs_running << "\n\n\n\n";
			GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			GlobalV::ofs_running << " |                                                                    |" << std::endl;
			GlobalV::ofs_running << " | Reading atom information in unitcell:                              |" << std::endl;
			GlobalV::ofs_running << " | From the input file and the structure file we know the number of   |" << std::endl;
			GlobalV::ofs_running << " | different elments in this unitcell, then we list the detail        |" << std::endl;
			GlobalV::ofs_running << " | information for each element, especially the zeta and polar atomic |" << std::endl;
			GlobalV::ofs_running << " | orbital number for each element. The total atom number is counted. |" << std::endl;
			GlobalV::ofs_running << " | We calculate the nearest atom distance for each atom and show the  |" << std::endl;
			GlobalV::ofs_running << " | Cartesian and Direct coordinates for each atom. We list the file   |" << std::endl;
			GlobalV::ofs_running << " | address for atomic orbitals. The volume and the lattice vectors    |" << std::endl;
			GlobalV::ofs_running << " | in real and reciprocal space is also shown.                        |" << std::endl;
			GlobalV::ofs_running << " |                                                                    |" << std::endl;
			GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			GlobalV::ofs_running << "\n\n\n\n";

			GlobalV::ofs_running << " READING UNITCELL INFORMATION" << std::endl;
			//========================
			// call read_atom_species
			//========================
			this->read_atom_species(ifa, GlobalV::ofs_running);

			//==========================
			// call read_atom_positions
			//==========================
			ok2 = this->read_atom_positions(ifa, GlobalV::ofs_running, GlobalV::ofs_warning);

			if(ok2)
			{
				for(int i=0;i<this->ntype;i++)
				{
					ModuleBase::Global_File::make_dir_atom( this->atoms[i].label );
				}
			}
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(ok);
	Parallel_Common::bcast_bool(ok2);
	if(GlobalV::NSPIN==4) 
	{
		Parallel_Common::bcast_bool(GlobalV::DOMAG);
		Parallel_Common::bcast_bool(GlobalV::DOMAG_Z);
	}
#endif
	if(!ok)
	{
		ModuleBase::WARNING_QUIT("UnitCell_pseudo::setup_cell","Can not find the file containing atom positions.!");
	}
	if(!ok2)
	{
		ModuleBase::WARNING_QUIT("UnitCell_pseudo::setup_cell","Something wrong during read_atom_positions.");
	}

#ifdef __MPI
	this->bcast_unitcell_pseudo();

	// mohan add 2010-09-29
	#ifdef __LCAO
	GlobalC::ORB.bcast_files(ntype, GlobalV::MY_RANK);
	#endif
#endif
	
	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, latvec must be read in.
	//========================================================
	assert(lat0 > 0.0);
	this->omega = abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
	if(this->omega<=0)
	{
		ModuleBase::WARNING_QUIT("setup_cell","omega <= 0 .");
	}
	else
	{
		GlobalV::ofs_running << std::endl;
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Volume (Bohr^3)", this->omega);
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Volume (A^3)", this->omega * pow(BOHR_TO_A, 3));
	}
		
	//==========================================================
	// Calculate recip. lattice vectors and dot products
	// latvec have the unit of lat0, but G has the unit 2Pi/lat0
	//==========================================================
	this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
	this->invGGT = GGT.Inverse();

    //LiuXh add 20180515
    this->GT0 = latvec.Inverse();
    this->G0  = GT.Transpose();
    this->GGT0 = G * GT;
    this->invGGT0 = GGT.Inverse();

	GlobalV::ofs_running << std::endl;
	outp.printM3(GlobalV::ofs_running,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec); 
	outp.printM3(GlobalV::ofs_running,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);
//	OUT(GlobalV::ofs_running,"lattice center x",latcenter.x);
//	OUT(GlobalV::ofs_running,"lattice center y",latcenter.y);
//	OUT(GlobalV::ofs_running,"lattice center z",latcenter.z);

	// read in non-local pseudopotential and ouput the projectors.

	GlobalV::ofs_running << "\n\n\n\n";
	GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " | Reading pseudopotentials files:                                    |" << std::endl;
	GlobalV::ofs_running << " | The pseudopotential file is in UPF format. The 'NC' indicates that |" << std::endl;
	GlobalV::ofs_running << " | the type of pseudopotential is 'norm conserving'. Functional of    |" << std::endl;
	GlobalV::ofs_running << " | exchange and correlation is decided by 4 given parameters in UPF   |" << std::endl;
	GlobalV::ofs_running << " | file.  We also read in the 'core correction' if there exists.      |" << std::endl;
	GlobalV::ofs_running << " | Also we can read the valence electrons number and the maximal      |" << std::endl;
	GlobalV::ofs_running << " | angular momentum used in this pseudopotential. We also read in the |" << std::endl;
	GlobalV::ofs_running << " | trail wave function, trail atomic density and local-pseudopotential|" << std::endl;
	GlobalV::ofs_running << " | on logrithmic grid. The non-local pseudopotential projector is also|" << std::endl;
	GlobalV::ofs_running << " | read in if there is any.                                           |" << std::endl;
	GlobalV::ofs_running << " |                                                                    |" << std::endl;
	GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	GlobalV::ofs_running << "\n\n\n\n";


	this->read_cell_pseudopots(s_pseudopot_dir);
	
	if(GlobalV::MY_RANK == 0) 
	{
		for(int it=0; it<this->ntype; it++)
		{
			Atom* atom = &atoms[it];
			std::stringstream ss;
			ss << GlobalV::global_out_dir << atom->label 
				<< "/" << atom->label
				<< ".NONLOCAL";
			std::ofstream ofs(ss.str().c_str());

			ofs << "<HEADER>" << std::endl;
			ofs << std::setw(10) << atom->label << "\t" << "label" << std::endl;
			ofs << std::setw(10) << atom->pp_type << "\t" << "pseudopotential type" << std::endl;
			ofs << std::setw(10) << atom->lmax << "\t" << "lmax" << std::endl;
			ofs << "</HEADER>" << std::endl;

			ofs << "\n<DIJ>" << std::endl;
			ofs << std::setw(10) << atom->nbeta << "\t" << "nummber of projectors." << std::endl;
			
			for(int ib=0; ib<atom->nbeta; ib++)
			{
				for(int ib2=0; ib2<atom->nbeta; ib2++)
				{
					ofs<<std::setw(10) << atom->lll[ib] 
						<< " " << atom->lll[ib2]
						<< " " << atom->dion(ib,ib2)<<std::endl;
				}
			}
			ofs << "</DIJ>" << std::endl;
			
			for(int i=0; i<atom->nbeta; i++)
			{
				ofs << "<PP_BETA>" << std::endl;
				ofs << std::setw(10) << i << "\t" << "the index of projectors." <<std::endl;
				ofs << std::setw(10) << atom->lll[i] << "\t" << "the angular momentum." <<std::endl;

				// mohan add
				// only keep the nonzero part.
				int cut_mesh = atom->mesh; 
				for(int j=atom->mesh-1; j>=0; --j)
				{
					if( abs( atom->betar(i,j) ) > 1.0e-10 )
					{
						cut_mesh = j; 
						break;
					}
				}
				if(cut_mesh %2 == 0) ++cut_mesh;

				ofs << std::setw(10) << cut_mesh << "\t" << "the number of mesh points." << std::endl;


				for(int j=0; j<cut_mesh; ++j)
				{
					ofs << std::setw(15) << atom->r[j]
						<< std::setw(15) << atom->betar(i, j)
						<< std::setw(15) << atom->rab[j] << std::endl;
				}
				ofs << "</PP_BETA>" << std::endl;
			}

			ofs.close();
		}
	}

#ifdef __MPI
	this->bcast_unitcell_pseudo2();
#endif	

	for(int it=0; it<ntype; it++)
	{
		for(int j=0; j<4; j++)
		{
			if(atoms[0].dft[j]!=atoms[it].dft[j])
			{
				GlobalV::ofs_warning << "\n type " << atoms[0].label << " functional is " 
				<< atoms[0].dft[0] << " " << atoms[0].dft[1] << " "
				<< atoms[0].dft[2] << " " << atoms[0].dft[3];
				
				GlobalV::ofs_warning << "\n type " << atoms[it].label << " functional is " 
				<< atoms[it].dft[0] << " " << atoms[it].dft[1] << " "
				<< atoms[it].dft[2] << " " << atoms[it].dft[3] << std::endl;
				
				ModuleBase::WARNING_QUIT("setup_cell","All DFT functional must consistent.");
			}
		}
	}

	// mohan add 2010-09-06
	// because the number of element type
	// will easily be ignored, so here
	// I warn the user again for each type.
	//for(int it=0; it<ntype; it++)
	//{
	//	GlobalC::xcf.which_dft(atoms[it].dft);
	//}

	// setup the total number of PAOs
	this->cal_natomwfc();

	// setup GlobalV::NLOCAL
	this->cal_nwfc();

	// setup GlobalV::NBANDS
	//this->cal_nelec();

	this->cal_meshx();

	// setup vdwd2 parameters
	//vdwd2_para.initset(*this);		// Peize Lin add 2021.03.09

//	std::stringstream ss;
//	ss << GlobalV::global_out_dir << "unitcell_pp.log";
//	print_unitcell_pseudo( ss.str() );
	return;
}

void UnitCell_pseudo::setup_cell_classic(
	const std::string &fn,
	std::ofstream &ofs_running,
	std::ofstream &ofs_warning)

{
	ModuleBase::TITLE("UnitCell_pseudo","setup_cell_classic");

	assert(ntype>0);

	// (1) init *Atom class array.
	this->atoms = new Atom[this->ntype];
	this->set_atom_flag = true;

	bool ok = true;
	bool ok2 = true;

	// (2) read in atom information
	if(GlobalV::MY_RANK == 0)
	{
		std::ifstream ifa(fn.c_str(), ios::in);
		if (!ifa)
		{
			ofs_warning << fn;
			ok = false;
		}

		if(ok)
		{
			ofs_running << "\n\n\n\n";
			ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
			ofs_running << " |                                                                    |" << std::endl;
			ofs_running << " | Reading atom information in unitcell for classic MD:               |" << std::endl;
			ofs_running << " | From the input file and the structure file we know the number of   |" << std::endl;
			ofs_running << " | different elments in this unitcell, then we list the detail        |" << std::endl;
			ofs_running << " | information for each element. The total atom number is counted.    |" << std::endl;
			ofs_running << " | We calculate the nearest atom distance for each atom and show the  |" << std::endl;
			ofs_running << " | Cartesian and Direct coordinates for each atom.                    |" << std::endl;
			ofs_running << " | The volume and the lattice vectors in real space is also shown.    |" << std::endl;
			ofs_running << " |                                                                    |" << std::endl;
			ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			ofs_running << "\n\n\n\n";

			ofs_running << " READING UNITCELL INFORMATION" << std::endl;
			//========================
			// call read_atom_species
			//========================
			this->read_atom_species(ifa, ofs_running);
			//==========================
			// call read_atom_positions
			//==========================
			ok2 = this->read_atom_positions(ifa, ofs_running, ofs_warning);
			if(ok2)
			{
				for(int i=0;i<this->ntype;i++)
				{
					ModuleBase::Global_File::make_dir_atom( this->atoms[i].label );
				}
			}
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(ok);
	Parallel_Common::bcast_bool(ok2);
#endif
	if(!ok)
	{
		ModuleBase::WARNING_QUIT("UnitCell_pseudo::setup_cell","Can not find the file containing atom positions.!");
	}
	if(!ok2)
	{
		ModuleBase::WARNING_QUIT("UnitCell_pseudo::setup_cell","Something wrong during read_atom_positions.");
	}

#ifdef __MPI
	this->bcast_unitcell();
#endif

	//========================================================
	// Calculate unit cell volume
	//========================================================
	assert(lat0 > 0.0);
	this->omega = abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
	if(this->omega<=0)
	{
		ModuleBase::WARNING_QUIT("setup_cell","omega <= 0 .");
	}
	else
	{
		ofs_running << std::endl;
		ModuleBase::GlobalFunc::OUT(ofs_running,"Volume (Bohr^3)", this->omega);
		ModuleBase::GlobalFunc::OUT(ofs_running,"Volume (A^3)", this->omega * pow(BOHR_TO_A, 3));
	}

	this->set_iat2itia();
}




//===========================================
// calculate the total number of local basis
// Target : nwfc, lmax,
// 			atoms[].stapos_wf
// 			GlobalV::NBANDS
//===========================================
void UnitCell_pseudo::cal_nwfc(void)
{
	ModuleBase::TITLE("UnitCell_pseudo","cal_nwfc");
	assert(ntype>0);
	assert(nat>0);

	//===========================
	// (1) set iw2l, iw2n, iw2m
	//===========================
	for(int it=0; it<ntype; it++)
	{
		this->atoms[it].set_index();
	}

	//===========================
	// (2) set namax and nwmax
	//===========================
	this->namax = 0;
	this->nwmax = 0;
	for(int it=0; it<ntype; it++)
	{
		this->namax = std::max( atoms[it].na, namax );
		this->nwmax = std::max( atoms[it].nw, nwmax );
	}
	assert(namax>0);
// for tests
//		OUT(GlobalV::ofs_running,"max input atom number",namax);
//		OUT(GlobalV::ofs_running,"max wave function number",nwmax);	

	//===========================
	// (3) set nwfc and stapos_wf
	//===========================
	GlobalV::NLOCAL = 0;
	for(int it=0; it<ntype; it++)
	{
		atoms[it].stapos_wf = GlobalV::NLOCAL;
		const int nlocal_it = atoms[it].nw * atoms[it].na;
		if(GlobalV::NSPIN!=4) 
		{
			GlobalV::NLOCAL += nlocal_it;
		}
		else 
		{
			GlobalV::NLOCAL += nlocal_it * 2;//zhengdy-soc
		}

// for tests
//		OUT(GlobalV::ofs_running,ss1.str(),nlocal_it);
//		OUT(GlobalV::ofs_running,"start position of local orbitals",atoms[it].stapos_wf);
	}
	
	//OUT(GlobalV::ofs_running,"NLOCAL",GlobalV::NLOCAL);
	GlobalV::ofs_running << " " << std::setw(40) << "NLOCAL" << " = " << GlobalV::NLOCAL <<std::endl;
	//========================================================
	// (4) set index for iat2it, iat2ia, itia2iat, itiaiw2iwt
	//========================================================

	this->set_iat2itia();
	
	//delete[] iat2ia;
	//this->iat2ia = new int[nat];// bug fix 2009-3-8

	// mohan add 2010-09-26
	assert(GlobalV::NLOCAL>0);
	delete[] iwt2iat;
	delete[] iwt2iw;
	this->iwt2iat = new int[GlobalV::NLOCAL];
	this->iwt2iw = new int[GlobalV::NLOCAL];

	this->itia2iat.create(ntype, namax);
	this->itiaiw2iwt.create(ntype, namax, nwmax*GlobalV::NPOL);
	int iat=0;
	int iwt=0;
	for(int it = 0;it < ntype;it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->itia2iat(it, ia) = iat;
			//this->iat2ia[iat] = ia;
			for(int iw=0; iw<atoms[it].nw * GlobalV::NPOL; iw++)
			{
				this->itiaiw2iwt(it, ia, iw) = iwt;
				this->iwt2iat[iwt] = iat;
				this->iwt2iw[iwt] = iw;
				++iwt;
			}
			++iat;
		}	
	}

	//========================
	// (5) set lmax and nmax
	//========================
	this->lmax = 0;	
	this->nmax = 0;
	for(int it=0; it<ntype; it++)
	{
		lmax = std::max( lmax, atoms[it].nwl);
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nmax = std::max( nmax, atoms[it].l_nchi[l] );
		}

		int nchi = 0;
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nchi += atoms[it].l_nchi[l];
		}
		this->nmax_total = std::max(nmax_total, nchi);
	}

	//=======================
	// (6) set lmax_ppwf
	//=======================
	this->lmax_ppwf = 0;
	for(int it=0; it< ntype; it++)
	{
		for(int ic=0; ic<atoms[it].nchi; ic++)
		{
			if( lmax_ppwf < atoms[it].lchi[ic] )
			{
				this->lmax_ppwf = atoms[it].lchi[ic]; 
			}
		}
	}

	/*
	for(int it=0; it< ntype; it++)
	{
		std::cout << " label=" << it << " nbeta=" << atoms[it].nbeta << std::endl;
		for(int ic=0; ic<atoms[it].nbeta; ic++)
		{
			std::cout << " " << atoms[it].lll[ic] << std::endl;
		}
	}
	*/

//	OUT(GlobalV::ofs_running,"lmax between L(pseudopotential)",lmax_ppwf);

	//=====================
	// Use localized basis
	//=====================
	if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
	{
		ModuleBase::GlobalFunc::AUTO_SET("NBANDS",GlobalV::NBANDS);
	}
	else // plane wave basis
	{
		//if(winput::after_iter && winput::sph_proj)
		//{
		//	if(GlobalV::NBANDS < GlobalV::NLOCAL)
		//	{
		//		ModuleBase::WARNING_QUIT("cal_nwfc","NBANDS must > GlobalV::NLOCAL !");
		//	}
		//}
	}

	return;
}

//======================
// Target : meshx
// Demand : atoms[].msh
//======================
void UnitCell_pseudo::cal_meshx()
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("UnitCell_pseudo","cal_meshx");
	this->meshx = 0;
	for (int it = 0;it < this->ntype;it++)
	{
		const int mesh = this->atoms[it].msh;
		if (mesh > this->meshx)
		{
			this->meshx = mesh;
		}
	}
	return;
}

//=========================
// Target : natomwfc
// Demand : atoms[].nchi
// 			atoms[].lchi
// 			atoms[].oc
// 			atoms[].na
//=========================
void UnitCell_pseudo::cal_natomwfc(void)
{
	if(GlobalV::test_pseudo_cell) ModuleBase::TITLE("UnitCell_pseudo","cal_natomwfc");

	this->natomwfc = 0;
	for (int it = 0;it < ntype;it++)
	{
		//============================
		// Use pseudo-atomic orbitals
		//============================
		int tmp=0;
		for (int l = 0;l < atoms[it].nchi;l++)
		{
			if (atoms[it].oc[l] >= 0)
			{
				if(GlobalV::NSPIN==4)
				{
					if(atoms[it].has_so)
					{
						tmp += 2 * atoms[it].lchi[l];
						if(fabs(atoms[it].jchi[l] - atoms[it].lchi[l] - 0.5)<1e-6)
						tmp += 2 ;
					}
					else
					{
						tmp += 2 * (2 * atoms[it].lchi[l] + 1);
					}
				}
				else
					tmp += 2 * atoms[it].lchi[l] + 1;
			}
		}
		natomwfc += tmp * atoms[it].na;
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"initial pseudo atomic orbital number",natomwfc);
	return;
}



//LiuXh add a new function here,
//20180515
void UnitCell_pseudo::setup_cell_after_vc(
        const std::string &s_pseudopot_dir,
		output &outp,
        const std::string &fn, std::ofstream &log)
{
    if(GlobalV::MY_RANK == 0)
    {
        //std::ifstream ifa(fn.c_str(), ios::in);
        //this->read_atom_species_after_vc(ifa);
    }

    assert(lat0 > 0.0);
    this->omega = abs(latvec.Det()) * this->lat0 * lat0 * lat0;
    if(this->omega <= 0)
    {
        ModuleBase::WARNING_QUIT("setup_cell_after_vc", "omega <= 0 .");
    }
    else
    {
        GlobalV::ofs_running << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Volume (Bohr^3)", this->omega);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Volume (A^3))", this->omega * pow(BOHR_TO_A, 3));
    }

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec has the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G  = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    for(int it=0; it<ntype; it++)
    {
        Atom* atom = &atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            atom->tau[ia] = atom->taud[ia] * latvec;
/*
#ifdef __MPI
Parallel_Common::bcast_double( atom->tau[ia].x );
Parallel_Common::bcast_double( atom->tau[ia].y );
Parallel_Common::bcast_double( atom->tau[ia].z );
Parallel_Common::bcast_double( atom->taud[ia].x );
Parallel_Common::bcast_double( atom->taud[ia].y );
Parallel_Common::bcast_double( atom->taud[ia].z );
#endif
*/
        }
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0;i<ntype;i++)
    {
        atoms[i].bcast_atom(); // bcast tau array
    }
#endif

    GlobalV::ofs_running << std::endl;
    outp.printM3(GlobalV::ofs_running,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec);
    outp.printM3(GlobalV::ofs_running,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);

    return;
}

//check if any atom can be moved
bool UnitCell_pseudo::if_atoms_can_move()const
{
	for(int it=0; it<this->ntype; it++)
    {
        Atom* atom = &atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            if(atom->mbl[ia].x||atom->mbl[ia].y||atom->mbl[ia].z) return 1;
		}
	}
	return 0;
}

//check if lattice vector can be changed
bool UnitCell_pseudo::if_cell_can_change()const
{
	//need to be fixed next
	return 1;
}