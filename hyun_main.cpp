#include "num_methods.h"
#include "matplotlibcpp.h"
#include <chrono>

namespace plt = matplotlibcpp;

void vec_sum(std::vector<float_T>& change, const std::vector<float_T>& add)
{
    size_t n = change.size();
    for (size_t i = 0; i < n; ++i) {
        change[i] += add[i];
    }
}

void vec_prod(std::vector<float_T>& change, const float_T& val)
{
    size_t n = change.size();
    for (size_t i = 0; i < n; ++i) {
        change[i] *= val;
    }
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t n;
    float_T D;
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();

    size_t N;
    std::cout << "Enter number of implementations\n";
    std::cin >> N;
    
    auto start_time = std::chrono::steady_clock::now();
    matrix sol, cum_sol = { std::vector<float_T>(n + 1), std::vector<float_T>(n + 1) };
    numeric_method solver(D);

    for (size_t i = 0; i < N; ++i) {
        sol = solver.stoch_hyun_method(n);
        vec_sum(cum_sol.x, sol.x);
    }
    cum_sol.t = sol.t;
    vec_prod(cum_sol.x, 1.0 / N);
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output << sol.t.size() << "\n";
    output.close();

    plt::named_plot("x(t)", cum_sol.t, cum_sol.x, "-");
    plt::grid(true);
    plt::show();

    return 0;
}