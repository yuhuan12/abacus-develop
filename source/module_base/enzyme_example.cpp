#include "matrix.h"
#include "math.h"
#include <vector>
#include <numeric>
#include <iostream>
#define D __enzyme_autodiff

using namespace std;

int enzyme_const;
int enzyme_out;
int enzyme_dup;

template<typename Return, typename ...T>
Return D(T...);

double wrapper_norm(matrix *m){ return m->norm(); }
double wrapper_max(matrix *m){ return m->max(); }
double wrapper_trace(matrix *m){ return m->trace_on(); }
double calc(double d){ return exp(d * d) / 2.0; };
double scale_sum(const vector<double> &vd, double scale) {
    return scale * accumulate(vd.begin(), vd.end(), 0.0);
}

int main() {
    cout << "Differentiating C standard library..." << endl;
    cout << D<double>(calc, 1.0) << endl << endl;
    cout << "Differentiating C++ standard library..." << endl;
    vector<double> vd(10, 1);
    double scale = 3;
    cout << D<double>(scale_sum,
        enzyme_const, vd,
        enzyme_out, scale
    ) << endl << endl;
    matrix m(2, 2);
    m(0, 0) = 1; m(0, 1) = 2; m(1, 0) = 3; m(1, 1) = 4;
    matrix _m(2, 2);
    cout << "Differentiating Fröbenius norm..." << endl;
    _m.zero_out();
    D<double>(wrapper_norm, &m, &_m);
    cout << _m << endl;
    cout << "Differentiating max..." << endl;
    _m.zero_out();
    D<double>(wrapper_max, &m, &_m);
    cout << _m << endl;
    cout << "Differentiating trace..." << endl;
    _m.zero_out();
    D<double>(wrapper_trace, &m, &_m);
    cout << _m << endl;
    return 0;
}
