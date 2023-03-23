#include <gsl/gsl_linalg.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <cmath>

typedef double float_T;

struct matrix
{
    std::vector<float_T> t;
    std::vector<float_T> x;
};

class numeric_method
{
private:
    float_T D;
public:
    numeric_method() { D = 1; }
    numeric_method(float_T D_) { D = D_; }

    void set_dispersion(float_T D_) { D = D_; }
    float_T get_dispersion() { return D; }

    matrix euler_method(float_T eps = pow(10, -6), float_T h = pow(10, -3));
    matrix stoch_euler_method(size_t N = 1000);
    matrix stoch_hyun_method(size_t N = 1000);
};
