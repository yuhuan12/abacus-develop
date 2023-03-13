#include "../global_file.h"
#include "../global_variable.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <fstream>
#include <unistd.h>

#ifdef __MPI
#include "mpi.h"
#endif

/************************************************
 *  unit test of functions in global_file.cpp
 ***********************************************/

/**
 * - Tested Function
 * - mkdiratom
 *   - generate atom dir for each type of atom
 * - openlog
 *   - Open the out file with the name *.log
 */

class GlobalFile : public testing::Test
{

};

TEST_F(GlobalFile,mkdirout)
{
	    std::string output;
	    testing::internal::CaptureStdout();
	    ModuleBase::Global_File::make_dir_out("Si","m",false,0,true,true);
	    output = testing::internal::GetCapturedStdout();
	    EXPECT_THAT(output,testing::HasSubstr("MAKE THE DIR"));
	    GlobalV::ofs_warning.close();
	    EXPECT_TRUE(GlobalV::ofs_running.is_open());
	    GlobalV::ofs_running.close();
	    std::string dd = "OUT.Si/running_m_1.log";
        remove(dd.c_str());

	    testing::internal::CaptureStdout();
	    ModuleBase::Global_File::make_dir_out("Si","md",false,0,true,false);
	    output = testing::internal::GetCapturedStdout();
	    EXPECT_THAT(output,testing::HasSubstr("MAKE THE STRU DIR"));
	    EXPECT_TRUE(GlobalV::ofs_running.is_open());
	    GlobalV::ofs_running.close();
	    GlobalV::ofs_warning.close();
	    std::string bb = "OUT.Si/running_md.log";
        remove(bb.c_str());

	    testing::internal::CaptureStdout();
	    ModuleBase::Global_File::make_dir_out("Si","md",true,0,true,true);
	    output = testing::internal::GetCapturedStdout();
	    EXPECT_THAT(output,testing::HasSubstr("MAKE THE MATRIX DIR"));
	    EXPECT_TRUE(GlobalV::ofs_running.is_open());
	    GlobalV::ofs_running.close();
	    std::string cc = "OUT.Si/running_md_1.log";
        remove(cc.c_str());
        std::string aa = "OUT.Si/warning.log";
        remove(aa.c_str());
        rmdir(GlobalV::global_stru_dir.c_str());
        rmdir(GlobalV::global_matrix_dir.c_str());
        rmdir(GlobalV::global_out_dir.c_str());
}

TEST_F(GlobalFile,mkdiratom)
{
        GlobalV::global_out_dir = "./";
        ModuleBase::Global_File::make_dir_atom("Si");
        int a = access("./Si/",0);
        EXPECT_EQ(a , 0);
        std::string ss = "./Si/";
        rmdir(ss.c_str());
}

TEST_F(GlobalFile,openlog)
{
        std::ofstream ofs;
        ModuleBase::Global_File::open_log(ofs,"Si","md",true);
        EXPECT_TRUE(ofs.is_open());
        ofs.close();
        ModuleBase::Global_File::open_log(ofs,"Si","md",false);
        EXPECT_TRUE(ofs.is_open());
        ofs.close();
        std::string sss = "Si.log"; 
        remove(sss.c_str());
}
