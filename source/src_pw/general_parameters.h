#ifndef GENERAL_PARAMETERS_H
#define GENERAL_PARAMETERS_H

#include <string>
namespace ABACUS
{
    class General_parameters
    {
        public:
            General_parameters()
            {
                calculation = "scf";
                basis_type = "pw";
                symmetry = 0;
                force = 0;
                stress = 0;
                //dft_plus_u = ""
            };
            ~General_parameters(){};

			/// "scf" : self consistent calculation.
			/// "nscf" : non-self consistent calculation.
			/// "relax" : cell relaxations, only ions can move. 
            /// "cell-relax" : cell relaxations, both ions and lattice vectors can move.
            /// "md" : molecular dynamics
			std::string calculation;
			/// "pw" : PW;
			/// "lcao_in_pw" : LCAO in pw;
			/// "lcao" :LCAO 
            std::string basis_type;
			bool symmetry;
            bool force;
            bool stress;
            bool dft_plus_u;
            std::string vdw_method;
            bool lspinorb;
            bool tddft;
            int nspin;
            double nelec;
            int nbands;
            int nbands_istate;
            std::string charge_extrap;
            std::string latname;
            int nbspline;
            std::string suffix;
            std::string dft_functional;
            double pseudo_rcut;
            bool mem_saver;
            double soc_lambda;
            std::string init_charge;
            std::string init_wfc;
            bool init_vel;
            bool init_mesh;
            int pw_seed;
            int k_par;
            int ntype;
    };

    class Input_path
    {
        public:
        Input_path()
        {
            kpoint_file = "";
            read_file_dir = "auto";
            stru_file = "";
        };
        ~Input_path(){};

        /// the name of file containing k points
        std::string kpoint_file;
        /// the filename of file containing atom positions
        std::string stru_file;
        /// directory of files for reading
        std::string read_file_dir;
    };
}

#endif
