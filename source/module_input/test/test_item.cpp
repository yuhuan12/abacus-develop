#include "../read_txt_input_list.h"
#include "../read_txt_input_process.h"
#include "gtest/gtest.h"
#include "../../abacus_parameters.h"

class item_test : public testing::Test
{
protected:
    

    void SetUp()
    {

        Read_Txt_Input::Input_List input_list;
        input_list.set_items();
        Read_Txt_Input::Input_Process input_process(input_list);
		input_process.read_and_convert("INPUT");
	}

	void test_input(std::string input_name)
	{
		Read_Txt_Input::Input_List input_list;
		input_list.set_items();
		Read_Txt_Input::Input_Process input_process(input_list);
		input_process.read_and_convert(input_name);
	}

    void TearDown()
    {
    }
};



TEST_F(item_test, mixing_beta)
{
    EXPECT_DOUBLE_EQ(ABACUS::para.scfp.mixing_beta, 0.5);
}
TEST_F(item_test, ecut)
{
	EXPECT_DOUBLE_EQ(ABACUS::para.pwp.ecut[0], 71.80);
}
TEST_F(item_test, scf_thr)
{
	EXPECT_DOUBLE_EQ(ABACUS::para.scfp.scf_thr, 3.5e-6);
}

TEST_F(item_test, check_output)
{
	test_input("INPUT-new-0");
	EXPECT_DOUBLE_EQ(ABACUS::para.scfp.mixing_beta, 0.5);
	EXPECT_DOUBLE_EQ(ABACUS::para.pwp.ecut[0], 71.80);
    EXPECT_DOUBLE_EQ(ABACUS::para.scfp.scf_thr, 3.5e-6);
}

int main(int argc, char **argv)
{
#ifdef __MPI
	MPI_Init(&argc, &argv);
#endif
	testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
	MPI_Finalize();
#endif
	return result;
}