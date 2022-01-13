#ifndef LCAO_PARAMETERS_H
#define LCAO_PARAMETERS_H

#include <string>
namespace ABACUS
{
class LCAO_parameters
{
public:
    LCAO_parameters()
    {
        nb2d = 1;
        lmaxmax = 2;
        lcao_dk = 0.01;
        lcao_dr = 0.01;
        lcao_rmax = 30;
        gamma_only = 0;
        bx = 2;
        by = 2;
        bz = 2;
        lcao_ecut = 0.0;
    };
    ~LCAO_parameters(){};
    int nb2d;
    int lmaxmax;
    double lcao_dk;
    double lcao_dr;
    double lcao_rmax;
    double lcao_ecut;
    bool gamma_only;
    ///
    /// big mesh ball
    ///
    int bx;
    int by;
    int bz;

};
}

#endif