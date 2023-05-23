#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <cmath>

#include "boost/random.hpp"

typedef double float_T;
typedef float_T* vec_T;

struct matrix
{
    vec_T t;
    vec_T x;
};

class numeric_method
{
private:
    float_T D; // Noise intensity
public:
    numeric_method() { D = 1; }
    numeric_method(float_T D_) { D = D_; }

    void set_dispersion(float_T D_) { D = D_; }
    float_T get_dispersion() { return D; }

    void euler_method(matrix& sol, size_t N, float_T a, float_T x0, 
        float_T tmin, float_T tmax, float_T RelTol = 1e-3, float_T AbsTol = 1e-6);
    void euler_maruyama_method(matrix& sol, size_t N, float_T a, float_T x0, 
        float_T tmin, float_T tmax);
    void hyun_method(matrix& sol, size_t N, float_T a, float_T x0, 
        float_T tmin, float_T tmax);
    void stoch_rk4_method(matrix& sol, size_t N, float_T a, float_T x0, 
        float_T tmin, float_T tmax);
};
