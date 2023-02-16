#include "gtest/gtest.h"
#include "gmock/gmock.h"

/************************************************
 *  unit test of read_pp
 ***********************************************/

/**
 * - Tested Functions:
 *   - ReadUPF100
 *     - read_pseudo_upf
 *       - read 1.0.0 version of upf pseudopotential file
 *     - read_pseudo_header
 *     - read_pseudo_mesh
 *     - read_pseudo_nl
 *     - read_pseudo_local
 *     - read_pseudo_pswfc
 *     - read_pseudo_rhoatom
 *     - read_pseudo_so
 *     - print_pseudo_upf
 *   - ReadUPF201
 *     - read_pseudo_upf201
 *       - read 2.0.1 version of upf pseudopotential file
 *     - getnameval
 *     - read_pseudo_upf201_r
 *     - read_pseudo_upf201_rab
 *     - read_pseudo_upf201_dij
 *     - read_pseudo_upf201_rhoatom
 *   - ReadUSppErr100
 *     - read_pseudo_nl
 *     - read_pseudo_nlcc
 *     - ultrasoft is not supported
 *   - HeaderErr201
 *     - ultrasoft and paw pp are not supported
 *   - ReadUPF201MESH2
 *     - a different "<PP_MESH" header in UPF file
 *   - ReadUPF201FR
 *     - read a full-relativistic pp
 *   - XCWarning
 *     - DFT functianal warning when mismatching between pp_file
 *     - and GlobalV::dft_functional happens
 */

#define private public
#include "module_cell/read_pp.h"

class ReadPPTest : public testing::Test
{
protected:
	std::string output;
};

TEST_F(ReadPPTest, ReadUPF100)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	ifs.open("./support/Te.pbe-rrkj.UPF");
	upf->read_pseudo_upf(ifs);
	EXPECT_TRUE(upf->has_so); // has soc info
	EXPECT_EQ(upf->nv,0); // number of version
	EXPECT_EQ(upf->psd,"Te"); // element label
	EXPECT_EQ(upf->pp_type,"NC"); // pp_type
	EXPECT_FALSE(upf->tvanp); // not ultrasoft
	EXPECT_FALSE(upf->nlcc); // no Nonlinear core correction
	EXPECT_EQ(upf->xc_func,"PBE"); // Exchange-Correlation functional
	EXPECT_EQ(upf->zp,6); // Z valence
	EXPECT_DOUBLE_EQ(upf->etotps,-15.54533017755); // total energy
	EXPECT_DOUBLE_EQ(upf->ecutwfc,0.0); // suggested cutoff for wfc
	EXPECT_DOUBLE_EQ(upf->ecutrho,0.0); // suggested cutoff for rho
	EXPECT_EQ(upf->lmax,2); // max angular momentum component
	EXPECT_EQ(upf->mesh,1191); // Number of points in mesh
	EXPECT_EQ(upf->nwfc,3); // Number of wavefunctions
	EXPECT_EQ(upf->nbeta,3); // Number of projectors
	EXPECT_EQ(upf->els[0],"5S"); // label for i-th atomic orbital
	EXPECT_EQ(upf->els[1],"5P"); // label for i-th atomic orbital
	EXPECT_EQ(upf->els[2],"5D"); // label for i-th atomic orbital
	EXPECT_EQ(upf->lchi[0],0); // angluar momentum of each atomic orbital
	EXPECT_EQ(upf->lchi[1],1); // angluar momentum of each atomic orbital
	EXPECT_EQ(upf->lchi[2],2); // angluar momentum of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->oc[0],2.0); // occupation of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->oc[1],3.0); // occupation of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->oc[2],1.0); // occupation of each atomic orbital
	EXPECT_DOUBLE_EQ(upf->r[0],1.75361916453E-05); // r
	EXPECT_DOUBLE_EQ(upf->r[upf->mesh-1],5.05901190442E+01); // r
	EXPECT_DOUBLE_EQ(upf->rab[0],2.19202395566E-07); // rab
	EXPECT_DOUBLE_EQ(upf->rab[upf->mesh-1],6.32376488053E-01); // rab
	EXPECT_DOUBLE_EQ(upf->vloc[0],-5.00890143222E+00); // vloc
	EXPECT_DOUBLE_EQ(upf->vloc[upf->mesh-1],-2.37200471955E-01); // vloc
	EXPECT_EQ(upf->lll[0],0); // BETA
	EXPECT_EQ(upf->kkbeta[0],957);
	EXPECT_DOUBLE_EQ(upf->beta(0,0),-1.82560984478E-03);
	EXPECT_DOUBLE_EQ(upf->beta(0,upf->kkbeta[0]-1),-1.61398366674E-03);
	EXPECT_EQ(upf->lll[1],0); // BETA
	EXPECT_EQ(upf->kkbeta[1],957);
	EXPECT_DOUBLE_EQ(upf->beta(1,0),1.51692136792E-03);
	EXPECT_DOUBLE_EQ(upf->beta(1,upf->kkbeta[1]-1),1.51458412218E-03);
	EXPECT_EQ(upf->lll[2],2); // BETA
	EXPECT_EQ(upf->kkbeta[2],957);
	EXPECT_DOUBLE_EQ(upf->beta(2,0),-3.10582746893E-13);
	EXPECT_DOUBLE_EQ(upf->beta(2,upf->kkbeta[2]-1),-4.17131335030E-04);
	EXPECT_EQ(upf->nd,4); // DIJ
	EXPECT_DOUBLE_EQ(upf->dion(0,0),-1.70394647943E-01);
	EXPECT_DOUBLE_EQ(upf->dion(0,1),-1.76521654672E-01);
	EXPECT_DOUBLE_EQ(upf->dion(1,1),-1.80323263809E-01);
	EXPECT_DOUBLE_EQ(upf->dion(2,2),1.16612440320E-02);
	EXPECT_DOUBLE_EQ(upf->chi(0,0),1.40252610787E-06);  // PSWFC
	EXPECT_DOUBLE_EQ(upf->chi(0,upf->mesh-1),6.15962544650E-25);
	EXPECT_DOUBLE_EQ(upf->chi(1,0),7.65306256201E-11);
	EXPECT_DOUBLE_EQ(upf->chi(1,upf->mesh-1),1.44320361049E-17);
	EXPECT_DOUBLE_EQ(upf->chi(2,0),4.37015997370E-16);
	EXPECT_DOUBLE_EQ(upf->chi(2,upf->mesh-1),6.20093850585E-05);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],3.93415898407E-12); // RhoAtom
	EXPECT_DOUBLE_EQ(upf->rho_at[upf->mesh-1],3.84516383534E-09);
	EXPECT_EQ(upf->nn[0],1); // nn
	EXPECT_EQ(upf->nn[1],2);
	EXPECT_EQ(upf->nn[2],3);
	EXPECT_DOUBLE_EQ(upf->jchi[0],0.0); // jchi
	EXPECT_DOUBLE_EQ(upf->jchi[1],0.0);
	EXPECT_NEAR(upf->jchi[3],0.0,1e-20);
	EXPECT_DOUBLE_EQ(upf->jjj[0],0.0); // jjj
	EXPECT_DOUBLE_EQ(upf->jjj[1],0.0);
	EXPECT_DOUBLE_EQ(upf->jjj[2],0.0);
	//EXPECT_EQ
	ifs.close();
	std::ofstream ofs;
	ofs.open("tmp");
	upf->print_pseudo_upf(ofs);
	ofs.close();
	ifs.open("tmp");
	getline(ifs, output);
	EXPECT_THAT(output, testing::HasSubstr("==== read_pseudo_upf ==="));
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, ReadUSppErr100)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	ifs.open("./support/Zn.pw91-n-van.UPF");
	//upf->read_pseudo_upf(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	// test output on screening
	// EXPECT_THAT(output,testing::HasSubstr("this function is called")); // read_pseudo_nlcc
	EXPECT_THAT(output,testing::HasSubstr("Ultra Soft Pseudopotential not available yet."));
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, ReadUPF201)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	ifs.open("./support/Cu_ONCV_PBE-1.0.upf");
	upf->read_pseudo_upf201(ifs);
	EXPECT_EQ(upf->psd,"Cu");
	EXPECT_EQ(upf->pp_type,"NC");
	EXPECT_FALSE(upf->has_so);
	EXPECT_FALSE(upf->nlcc);
	EXPECT_EQ(upf->xc_func,"PBE");
	EXPECT_EQ(upf->zp,19);
	EXPECT_DOUBLE_EQ(upf->etotps,-1.82394100797E+02);
	EXPECT_EQ(upf->lmax,2);
	EXPECT_EQ(upf->mesh,601); // mesh -= 1 at line 388 (why? Let's see)
	EXPECT_EQ(upf->nwfc,0);
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_DOUBLE_EQ(upf->r[0],0.0);
	EXPECT_DOUBLE_EQ(upf->r[601],6.01);
	EXPECT_DOUBLE_EQ(upf->rab[0],0.01);
	EXPECT_DOUBLE_EQ(upf->rab[601],0.01);
	EXPECT_DOUBLE_EQ(upf->vloc[0],-5.3426582174E+01);
	EXPECT_DOUBLE_EQ(upf->vloc[601],-6.3227960060E+00);
	EXPECT_EQ(upf->lll[0],0);
	EXPECT_EQ(upf->kkbeta[0],196);
	EXPECT_DOUBLE_EQ(upf->beta(0,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(0,601),0.0);
	EXPECT_EQ(upf->lll[1],0);
	EXPECT_EQ(upf->kkbeta[1],196);
	EXPECT_DOUBLE_EQ(upf->beta(1,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(1,601),0.0);
	EXPECT_EQ(upf->lll[2],1);
	EXPECT_EQ(upf->kkbeta[2],196);
	EXPECT_DOUBLE_EQ(upf->beta(2,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(2,601),0.0);
	EXPECT_EQ(upf->lll[3],1);
	EXPECT_EQ(upf->kkbeta[3],196);
	EXPECT_DOUBLE_EQ(upf->beta(3,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(3,601),0.0);
	EXPECT_EQ(upf->lll[4],2);
	EXPECT_EQ(upf->kkbeta[4],196);
	EXPECT_DOUBLE_EQ(upf->beta(4,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(4,601),0.0);
	EXPECT_EQ(upf->lll[5],2);
	EXPECT_EQ(upf->kkbeta[5],196);
	EXPECT_DOUBLE_EQ(upf->beta(5,0),0.0);
	EXPECT_DOUBLE_EQ(upf->beta(5,601),0.0);
	EXPECT_DOUBLE_EQ(upf->dion(0,0),-6.6178420255E+00);
	EXPECT_DOUBLE_EQ(upf->dion(5,5),-7.0938557228E+00);
	EXPECT_DOUBLE_EQ(upf->rho_at[0],0.0);
	EXPECT_DOUBLE_EQ(upf->rho_at[601],3.1742307110E-02);
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, HeaderErr2011)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// 1st
	ifs.open("./support/HeaderError1");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf201(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("unknown pseudo type"));
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, HeaderErr2012)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// 2nd
	ifs.open("./support/HeaderError2");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf201(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("ULTRASOFT PSEUDOPOTENTIAL IS NOT SUPPORTED"));
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, HeaderErr2013)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// 3rd
	ifs.open("./support/HeaderError3");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_EXIT(upf->read_pseudo_upf201(ifs),
			::testing::ExitedWithCode(0),"");
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("PAW PSEUDOPOTENTIAL IS NOT SUPPORTED"));
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, HeaderErr2014)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// 4th
	GlobalV::ofs_warning.open("warning.log");
	ifs.open("./support/HeaderError4");
	upf->read_pseudo_upf201(ifs);
	GlobalV::ofs_warning.close();
	ifs.close();
	ifs.open("warning.log");
	getline(ifs,output);
	EXPECT_THAT(output,testing::HasSubstr("arbitrary is not read in. Please add this parameter in read_pp_upf201.cpp if needed."));
	ifs.close();
	remove("warning.log");
	delete upf;
}

TEST_F(ReadPPTest, ReadUPF201FR)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// this is a dojo full-relativisitic pp
	ifs.open("./support/C.upf");
	upf->read_pseudo_upf201(ifs);
	EXPECT_EQ(upf->psd,"C");
	EXPECT_TRUE(upf->has_so);
	EXPECT_TRUE(upf->nlcc);
	EXPECT_EQ(upf->mesh,1247);
	//RELBETA
	EXPECT_EQ(upf->nbeta,6);
	EXPECT_EQ(upf->lll[0],0);
	EXPECT_EQ(upf->lll[1],0);
	EXPECT_EQ(upf->lll[2],1);
	EXPECT_EQ(upf->lll[3],1);
	EXPECT_EQ(upf->lll[4],1);
	EXPECT_EQ(upf->lll[5],1);
	EXPECT_DOUBLE_EQ(upf->jjj[0],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[1],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[2],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[3],1.5);
	EXPECT_DOUBLE_EQ(upf->jjj[4],0.5);
	EXPECT_DOUBLE_EQ(upf->jjj[5],1.5);
	//RELWFC
	EXPECT_EQ(upf->nwfc,3);
	EXPECT_EQ(upf->nn[0],1);
	EXPECT_EQ(upf->nn[1],2);
	EXPECT_EQ(upf->nn[2],2);
	EXPECT_EQ(upf->lchi[0],0);
	EXPECT_EQ(upf->lchi[1],1);
	EXPECT_EQ(upf->lchi[2],1);
	EXPECT_DOUBLE_EQ(upf->jchi[0],0.5);
	EXPECT_DOUBLE_EQ(upf->jchi[1],1.5);
	EXPECT_DOUBLE_EQ(upf->jchi[2],0.5);
	//PSWFC
	EXPECT_EQ(upf->els[0],"2S");
	EXPECT_EQ(upf->lchi[0],0);
	EXPECT_DOUBLE_EQ(upf->oc[0],2.0);
	EXPECT_EQ(upf->els[1],"2P");
	EXPECT_EQ(upf->lchi[1],1);
	EXPECT_DOUBLE_EQ(upf->oc[1],1.333);
	EXPECT_EQ(upf->els[2],"2P");
	EXPECT_EQ(upf->lchi[2],1);
	EXPECT_DOUBLE_EQ(upf->oc[2],0.667);
	EXPECT_DOUBLE_EQ(upf->chi(0,0),2.0715339166E-12);
	EXPECT_DOUBLE_EQ(upf->chi(2,upf->mesh-1),1.1201306967E-03);
	//NLCC
	EXPECT_DOUBLE_EQ(upf->rho_atc[0],8.7234550809E-01);
	EXPECT_DOUBLE_EQ(upf->rho_atc[upf->mesh-1],0.0);
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, ReadUPF201MESH2)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// this pp file has gipaw, thus a different header
	ifs.open("./support/Fe.pbe-sp-mt_gipaw.UPF");
	upf->read_pseudo_upf201(ifs);
	EXPECT_EQ(upf->psd,"Fe");
	ifs.close();
	delete upf;
}

TEST_F(ReadPPTest, XCWarning)
{
	Pseudopot_upf* upf = new Pseudopot_upf;
	std::ifstream ifs;
	// this pp file has gipaw, thus a different header
	// dft_functional warning
	GlobalV::DFT_FUNCTIONAL="LDA";
	ifs.open("./support/Fe.pbe-sp-mt_gipaw.UPF");
	//upf->read_pseudo_upf201(ifs);
	testing::internal::CaptureStdout();
	EXPECT_NO_THROW(upf->read_pseudo_upf201(ifs));
	output = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output,testing::HasSubstr("dft_functional readin is: LDA"));
	EXPECT_THAT(output,testing::HasSubstr("dft_functional in pseudopot file is: PBE"));
	ifs.close();
	delete upf;
}
#undef private
