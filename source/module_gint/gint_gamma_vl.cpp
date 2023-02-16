//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "module_orbital/ORB_read.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/blas_connector.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

#include "module_hamilt_lcao/hamilt_lcaodft/global_fp.h" // mohan add 2021-01-30

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
}

void Gint_Gamma::cal_vlocal(Gint_inout *inout, const bool new_e_iteration)
{
	const int max_size = GlobalC::GridT.max_atom;
	const int lgd = GlobalC::GridT.lgd;

	if(inout->job==Gint_Tools::job_type::vlocal || inout->job==Gint_Tools::job_type::vlocal_meta)
	{
        if (max_size >0 && lgd > 0)
        {
            pvpR_grid = new double[lgd*lgd];
            ModuleBase::GlobalFunc::ZEROS(pvpR_grid, lgd*lgd);            
        }

        this->cal_gint(inout);

		this->vl_grid_to_2D(lgd,inout->lm[0], new_e_iteration);

		if (max_size >0 && lgd > 0) delete[] pvpR_grid;
	}
}

#ifdef __MPI
//------------------------------------------------------------------
// mohan add notes: 2021-03-11
// this subroutine is used to transform data from grid integrals
// to 2D-block distribution
// s stands for 'sender' and r stands for 'receiver'
//------------------------------------------------------------------
inline int setBufferParameter(
	MPI_Comm comm_2D,
	int blacs_ctxt,
	int nblk,
	int& s_index_siz,
	int*& s_local_index,
	int*& s_siz_pro,
	int*& s_dis_pro,
	int& s_siz,
	double*& s_buffer,
	int& r_index_siz,
	int*& r_global_index,
	int*& r_siz_pro,
	int*& r_dis_pro,
	int& r_siz,
	double*& r_buffer)
{
	//-----------------------------------------
    // setup blacs parameters
	//-----------------------------------------
    int nprows, npcols, nprocs;
    int myprow, mypcol, myproc;

    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

	//-----------------------------------------
	// set index of current proor: myproc
	// set number of total proors: nprocs
	//-----------------------------------------
    Cblacs_pinfo(&myproc, &nprocs);

    // initialize data arrays
    delete[] s_siz_pro;
    delete[] s_dis_pro;
    delete[] r_siz_pro;
    delete[] r_dis_pro;

    s_siz_pro=new int[nprocs];
    s_dis_pro=new int[nprocs];
    r_siz_pro=new int[nprocs];
    r_dis_pro=new int[nprocs];

	//---------------------------------------------------------------------
    // build the local index to be sent to other pro (s_local_index),
    // the global index to be received from other pro (r_global_index),
    // the send/receive siz/dis for data exchange by MPI_Alltoall
	//---------------------------------------------------------------------
    s_index_siz=GlobalC::GridT.lgd*GlobalC::GridT.lgd*2;

    delete[] s_local_index;
    s_local_index=new int[s_index_siz];

    int *s_global_index=new int[s_index_siz];

    int pos=0;
    s_siz_pro[0]=0;
    for(int iproc=0; iproc<nprocs; ++iproc)
    {
        s_dis_pro[iproc]=pos;

        int iprow=0;
		int ipcol=0;
        Cblacs_pcoord(blacs_ctxt, iproc, &iprow, &ipcol);

        // find out the global index and local index of elements
		// in each pro based on 2D block cyclic distribution
        for(int irow=0, grow=0; grow<GlobalV::NLOCAL; ++irow)
        {
            grow=Local_Orbital_wfc::globalIndex(irow, nblk, nprows, iprow);
            if (grow >= GlobalV::NLOCAL)
                continue;
            int lrow = GlobalC::GridT.trace_lo[grow];
            if (lrow < 0)
                continue;

            for(int icol=0, gcol=0; gcol<GlobalV::NLOCAL; ++icol)
            {
                gcol=Local_Orbital_wfc::globalIndex(icol,nblk, npcols, ipcol);
                if (gcol >= GlobalV::NLOCAL)
                    continue;
                int lcol = GlobalC::GridT.trace_lo[gcol];
                if (lcol < 0)
                    continue;
                s_global_index[pos]=grow;
                s_global_index[pos+1]=gcol;
                s_local_index[pos]=lrow;
                s_local_index[pos+1]=lcol;
                pos+=2;
            }
        }
        s_siz_pro[iproc]=pos-s_dis_pro[iproc];
    }

    MPI_Alltoall(s_siz_pro, 1, MPI_INT,
                 r_siz_pro, 1, MPI_INT, comm_2D);

    r_index_siz=r_siz_pro[0];
    r_dis_pro[0]=0;
    for(int i=1; i<nprocs; ++i)
    {
        r_index_siz+=r_siz_pro[i];
        r_dis_pro[i]=r_dis_pro[i-1]+r_siz_pro[i-1];
    }

	delete[] r_global_index;
	r_global_index=new int[r_index_siz];

    // send the global index in sendBuffer to recvBuffer
    MPI_Alltoallv(s_global_index, s_siz_pro, s_dis_pro, MPI_INT,
                  r_global_index, r_siz_pro, r_dis_pro, MPI_INT, comm_2D);

    delete [] s_global_index;

    // the s_siz_pro, s_dis_pro, r_siz_pro,
    // and r_dis_pro will be used in transfer s_buffer, which
    // is half siz of s_global_index
    // we have to rebuild the siz and dis for each pro
    for (int iproc=0; iproc < nprocs; ++iproc)
    {
        s_siz_pro[iproc]=s_siz_pro[iproc]/2;
        s_dis_pro[iproc]=s_dis_pro[iproc]/2;
        r_siz_pro[iproc]=r_siz_pro[iproc]/2;
        r_dis_pro[iproc]=r_dis_pro[iproc]/2;
    }

    s_siz=s_index_siz/2;
	delete[] s_buffer;
	s_buffer=new double[s_siz];

    r_siz=r_index_siz/2;
	delete[] r_buffer;
	r_buffer=new double[r_siz];

    return 0;
}
#endif

void Gint_Gamma::vl_grid_to_2D(const int lgd_now, LCAO_Matrix &lm, const bool new_e_iteration)
{
    ModuleBase::timer::tick("Gint_Gamma","distri_vl");
    // setup send buffer and receive buffer size
    // OUT(GlobalV::ofs_running, "Start transforming vlocal from grid distribute to 2D block");
    if(new_e_iteration)
    {
        ModuleBase::timer::tick("Gint_Gamma","distri_vl_index");
        #ifdef __MPI
        setBufferParameter(lm.ParaV->comm_2D, lm.ParaV->blacs_ctxt, lm.ParaV->nb,
                           this->sender_index_size, this->sender_local_index,
                           this->sender_size_process, this->sender_displacement_process,
                           this->sender_size, this->sender_buffer,
                           this->receiver_index_size, this->receiver_global_index,
                           this->receiver_size_process, this->receiver_displacement_process,
                           this->receiver_size, this->receiver_buffer);
        #endif
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "vlocal exchange index is built");
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "buffer size(M):", (this->sender_size+this->receiver_size)*sizeof(double)/1024/1024);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "buffer index size(M):", (this->sender_index_size+this->receiver_index_size)*sizeof(int)/1024/1024);
        ModuleBase::timer::tick("Gint_Gamma","distri_vl_index");
    }

    ModuleBase::timer::tick("Gint_Gamma","distri_vl_value");

    // put data to send buffer
    for(int i=0; i<this->sender_index_size; i+=2)
    {
        const int irow=this->sender_local_index[i];
        const int icol=this->sender_local_index[i+1];
        if(irow<=icol)
		{
            this->sender_buffer[i/2]=pvpR_grid[irow*lgd_now+icol];
		}
        else
		{
            this->sender_buffer[i/2]=pvpR_grid[icol*lgd_now+irow];
		}
    }
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "vlocal data are put in sender_buffer, size(M):", this->sender_size*8/1024/1024);

    // use mpi_alltoall to get local data
    #ifdef __MPI
    MPI_Alltoallv(this->sender_buffer, this->sender_size_process, this->sender_displacement_process, MPI_DOUBLE,
                  this->receiver_buffer, this->receiver_size_process,
					this->receiver_displacement_process, MPI_DOUBLE, lm.ParaV->comm_2D);
    #endif
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "vlocal data are exchanged, received size(M):", this->receiver_size*8/1024/1024);

    // put local data to H matrix
    for(int i=0; i<this->receiver_index_size; i+=2)
    {
        const int g_row=this->receiver_global_index[i];
        const int g_col=this->receiver_global_index[i+1];
        lm.set_HSgamma(g_row,g_col,this->receiver_buffer[i/2],'L', lm.Hloc.data());
    }

    ModuleBase::timer::tick("Gint_Gamma","distri_vl_value");
    ModuleBase::timer::tick("Gint_Gamma","distri_vl");
}
