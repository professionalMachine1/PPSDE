#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <cmath>

typedef double float_T;
typedef std::vector<float_T> vec_T;

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

    matrix euler_method(float_T RelTol = 1e-3, float_T AbsTol = 1e-6);
    void euler_maruyama_method(matrix& sol, size_t N);
    void hyun_method(matrix& sol, size_t N);
    void stoch_rk4_method(matrix& sol, size_t N);
};
