#include "matrix.h"
#include <iostream>
#define grad __enzyme_autodiff

void grad(double (*)(matrix *), matrix *, matrix *);

int main() {
    matrix m(2, 2);
    m(0, 0) = 1; m(0, 1) = 2; m(1, 0) = 3; m(1, 1) = 4;
    matrix _m(2, 2);
    std::cout << "Differentiating Fröbenius norm..." << std::endl;
    _m.zero_out();
    grad([](matrix *m){ return m->norm(); }, &m, &_m);
    std::cout << _m << std::endl;
    std::cout << "Differentiating max..." << std::endl;
    _m.zero_out();
    grad([](matrix *m){ return m->max(); }, &m, &_m);
    std::cout << _m << std::endl;
    std::cout << "Differentiating trace..." << std::endl;
    _m.zero_out();
    grad([](matrix *m){ return m->trace_on(); }, &m, &_m);
    std::cout << _m << std::endl;
    return 0;
}
