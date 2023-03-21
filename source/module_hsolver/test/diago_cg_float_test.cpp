#include "module_base/inverse_matrix.h"
#include "module_base/lapack_connector.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_psi/psi.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "../diago_cg.h"
#include "../diago_iter_assist.h"
#include "diago_mock.h"
#include "mpi.h"
#include "module_pw/test/test_tool.h"
#include <complex>

#include "gtest/gtest.h"
#include <random>

/************************************************
 *  unit test of functions in Diago_CG
 ***********************************************/

/**
 * Class Diago_CG is an approach for eigenvalue problems
 * This unittest test the function Diago_CG::diag() for FPTYPE=float, Device=cpu
 * with different examples.
 *  - the Hermite matrices (npw=500,1000) produced using random numbers and with sparsity of 0%, 60%, 80%
 *  - the Hamiltonian matrix read from "data-H", produced by using out_hs in INPUT of a LCAO calculation
 *  - a 2x2 Hermite matrix for learning and checking
 *
 * Note:
 * The test is passed when the eignvalues are closed to these calculated by LAPACK.
 * It is used together with a header file diago_mock.h.
 * The default Hermite matrix generated here is real symmetric, one can add an imaginary part
 * by changing two commented out lines in diago_mock.h.
 *
 */

// call lapack in order to compare to cg
void lapackEigen(int &npw, std::vector<std::complex<float>> &hm, float *e, bool outtime = false)
{
    clock_t start, end;
    start = clock();
    int lwork = 2 * npw;
    std::complex<float> *work2 = new std::complex<float>[lwork];
    float *rwork = new float[3 * npw - 2];
    int info = 0;
    char tmp_c1 = 'V', tmp_c2 = 'U';
    cheev_(&tmp_c1, &tmp_c2, &npw, hm.data(), &npw, e, work2, &lwork, rwork, &info);
    end = clock();
    if (outtime)
        std::cout << "Lapack Run time: " << (float)(end - start) / CLOCKS_PER_SEC << " S" << std::endl;
    delete[] rwork;
    delete[] work2;
}

class DiagoCGPrepare
{
  public:
    DiagoCGPrepare(int nband, int npw, int sparsity, bool reorder, float eps, int maxiter, float threshold)
        : nband(nband), npw(npw), sparsity(sparsity), reorder(reorder), eps(eps), maxiter(maxiter),
          threshold(threshold)
    {
#ifdef __MPI	
		MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
#endif	
    }

    int nband, npw, sparsity, maxiter, notconv;
    // eps is the convergence threshold within cg_diago
    float eps, avg_iter;
    bool reorder;
    float threshold;
    int nprocs=1, mypnum=0;
    // threshold is the comparison standard between cg and lapack

    void CompareEigen(float *precondition)
    {
        // calculate eigenvalues by LAPACK;
        float *e_lapack = new float[npw];
        auto ev = DIAGOTEST::hmatrix_f;

        if(mypnum == 0)  lapackEigen(npw, ev, e_lapack, false);

        // initial guess of psi by perturbing lapack psi
        std::vector<std::complex<float>> psiguess(nband * npw);
        std::default_random_engine p(1);
        std::uniform_int_distribution<unsigned> u(1, 10);

        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
		        float rand = static_cast<float>(u(p))/10.;
                // psiguess(i,j) = ev(j,i)*(1+rand);
                psiguess[i * npw + j] = ev[j * DIAGOTEST::h_nc + i] * rand;
            }
        }

        // run cg
	//======================================================================
        float *en = new float[npw];
        int ik = 1;
	    hamilt::Hamilt<float>* ha;
	    ha =new hamilt::HamiltPW<float>(nullptr);
	    int* ngk = new int [1];
	    //psi::Psi<std::complex<float>> psi(ngk,ik,nband,npw);
	    psi::Psi<std::complex<float>> psi;
	    psi.resize(ik,nband,npw);
	    //psi.fix_k(0);
        for (int i = 0; i < nband; i++)
        {
            for (int j = 0; j < npw; j++)
            {
	            psi(i,j)=psiguess[i * npw + j];
	        }
	    }	

        psi::Psi<std::complex<float>> psi_local;
        float* precondition_local;
        DIAGOTEST::npw_local = new int[nprocs];
#ifdef __MPI				
	    DIAGOTEST::cal_division(DIAGOTEST::npw);
	    DIAGOTEST::divide_hpsi_f(psi,psi_local); //will distribute psi and Hmatrix to each process
	    precondition_local = new float[DIAGOTEST::npw_local[mypnum]];
	    DIAGOTEST::divide_psi_f<float>(precondition,precondition_local);	
#else
	    DIAGOTEST::hmatrix_local_f = DIAGOTEST::hmatrix_f;
	    DIAGOTEST::npw_local[0] = DIAGOTEST::npw;
	    psi_local = psi;
	    precondition_local = new float[DIAGOTEST::npw];
	    for(int i=0;i<DIAGOTEST::npw;i++) precondition_local[i] = precondition[i];
#endif
        hsolver::DiagoCG<float> cg(precondition_local);
        psi_local.fix_k(0);
        float start, end;
        start = MPI_Wtime();
        cg.diag(ha,psi_local,en); 
        end = MPI_Wtime();
        //if(mypnum == 0) printf("diago time:%7.3f\n",end-start);
        delete [] DIAGOTEST::npw_local;
	    delete [] precondition_local;
	    //======================================================================
        for (int i = 0; i < nband; i++)
        {
            EXPECT_NEAR(en[i], e_lapack[i], threshold)<<i;
        }

        delete[] en;
        delete[] e_lapack;
        delete ha;
    }
};

class DiagoCGFloatTest : public ::testing::TestWithParam<DiagoCGPrepare>
{
};

TEST_P(DiagoCGFloatTest, RandomHamilt)
{
    DiagoCGPrepare dcp = GetParam();
    //std::cout << "npw=" << dcp.npw << ", nband=" << dcp.nband << ", sparsity="
    //		  << dcp.sparsity << ", eps=" << dcp.eps << std::endl;
    hsolver::DiagoIterAssist<float>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<float>::PW_DIAG_THR = dcp.eps;
    //std::cout<<"maxiter "<<hsolver::DiagoIterAssist<float>::PW_DIAG_NMAX<<std::endl;
    //std::cout<<"eps "<<hsolver::DiagoIterAssist<float>::PW_DIAG_THR<<std::endl;
    HPsi_f hpsi(dcp.nband, dcp.npw, dcp.sparsity);
    DIAGOTEST::hmatrix_f = hpsi.hamilt();

    DIAGOTEST::npw = dcp.npw;
    // ModuleBase::ComplexMatrix psi = hpsi.psi();
    dcp.CompareEigen(hpsi.precond());
}

INSTANTIATE_TEST_SUITE_P(VerifyCG,
                         DiagoCGFloatTest,
                         ::testing::Values(
                             // nband, npw, sparsity, reorder, eps, maxiter, threshold
                             DiagoCGPrepare(10, 200, 0, true, 1e-6, 300, 1e-0),
                             DiagoCGPrepare(10, 200, 6, true, 1e-6, 300, 1e-0),
                             DiagoCGPrepare(10, 400, 8, true, 1e-6, 300, 1e-0),
                             DiagoCGPrepare(10, 600, 8, true, 1e-6, 300, 1e-0))); 
                            //DiagoCGPrepare(40, 2000, 8, true, 1e-5, 500, 1e-2))); 
			    // the last one is passed but time-consumming.

// check that the mock class HPsi work well
// in generating a Hermite matrix
TEST(DiagoCGFloatTest, Hamilt)
{
    int dim = 2;
    int nbnd = 2;
    HPsi_f hpsi(nbnd, dim);
    std::vector<std::complex<float>> hm = hpsi.hamilt();
    EXPECT_EQ(DIAGOTEST::h_nr, 2);
    EXPECT_EQ(DIAGOTEST::h_nc, 2);
    EXPECT_EQ(hm[0].imag(), 0.0);
    EXPECT_EQ(hm[DIAGOTEST::h_nc + 1].imag(), 0.0);
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).real(), hm[1].real());
    EXPECT_EQ(conj(hm[DIAGOTEST::h_nc]).imag(), hm[1].imag());
}

// check that lapack work well
// for an eigenvalue problem
/*TEST(DiagoCGFloatTest, ZHEEV)
{
    int dim = 100;
    int nbnd = 2;
    HPsi hpsi(nbnd, dim);
    std::vector<std::complex<float>> hm = hpsi.hamilt();
    std::vector<std::complex<float>> hm_backup = hm;
    ModuleBase::ComplexMatrix eig(dim, dim);
    float e[dim];
    // using zheev to do a direct test
    lapackEigen(dim, hm, e);
    eig = transpose(hm, true) * hm_backup * hm;
    // for (int i=0;i<dim;i++) std::cout<< " e[i] "<<e[i]<<std::endl;
    for (int i = 0; i < dim; i++)
    {
        EXPECT_NEAR(e[i], eig(i, i).real(), 1e-10);
    }
}*/

// cg for a 2x2 matrix
#ifdef __MPI
#else
TEST(DiagoCGFloatTest, TwoByTwo)
{
    int dim = 2;
    int nband = 2;
    ModuleBase::ComplexMatrix hm(2, 2);
    hm(0, 0) = std::complex<float>{4.0, 0.0};
    hm(0, 1) = std::complex<float>{1.0, 0.0};
    hm(1, 0) = std::complex<float>{1.0, 0.0};
    hm(1, 1) = std::complex<float>{3.0, 0.0};
    // nband, npw, sub, sparsity, reorder, eps, maxiter, threshold
    DiagoCGPrepare dcp(nband, dim, 0, true, 1e-4, 50, 1e-0);
    hsolver::DiagoIterAssist<float>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<float>::PW_DIAG_THR = dcp.eps;
    HPsi hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix_f = hm;
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}
#endif

TEST(DiagoCGFloatTest, readH)
{
    // read Hamilt matrix from file data-H
    std::vector<std::complex<float>> hm;
    std::ifstream ifs;
    ifs.open("H-KPoints-Si64.dat");
    DIAGOTEST::readh(ifs, hm);
    ifs.close();
    int dim = DIAGOTEST::npw;
    int nband = 10; // not nband < dim, here dim = 26 in data-H
    // nband, npw, sparsity, reorder, eps, maxiter, threshold
    DiagoCGPrepare dcp(nband, dim, 0, true, 1e-5, 500, 1e-0);
    hsolver::DiagoIterAssist<float>::PW_DIAG_NMAX = dcp.maxiter;
    hsolver::DiagoIterAssist<float>::PW_DIAG_THR = dcp.eps;
    HPsi_f hpsi;
    hpsi.create(nband, dim);
    DIAGOTEST::hmatrix_f = hpsi.hamilt();
    DIAGOTEST::npw = dim;
    dcp.CompareEigen(hpsi.precond());
}

int main(int argc, char **argv)
{
	int nproc = 1, myrank = 0;

#ifdef __MPI
	int nproc_in_pool, kpar=1, mypool, rank_in_pool;
    setupmpi(argc,argv,nproc, myrank);
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
    GlobalV::NPROC_IN_POOL = nproc;
#else
	MPI_Init(&argc, &argv);	
#endif

    testing::InitGoogleTest(&argc, argv);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (myrank != 0) delete listeners.Release(listeners.default_result_printer());

    int result = RUN_ALL_TESTS();
    if (myrank == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
	}

    MPI_Finalize();
	return 0;
}
