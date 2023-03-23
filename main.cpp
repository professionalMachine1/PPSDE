#include "num_methods.h"
#include "matplotlibcpp.h"
#include <chrono>

namespace plt = matplotlibcpp;

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    std::cout << "Choose an integration method: \n";
    std::cout << " 1 - euler method\n 2 - stochastic euler method\n";
    std::cout << " 3 - stochastic hyun method\n";

    int32_t method_number;
    std::cin >> method_number;

    float_T n;
    float_T D;
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();
    
    auto start_time = std::chrono::steady_clock::now();
    matrix sol;
    numeric_method solver(D);

    switch (method_number)
    {
    case 1:
        sol = solver.euler_method();
        break;
    case 2:
        sol = solver.stoch_euler_method(n);
        break;
    case 3:
        sol = solver.stoch_hyun_method(n);
        break;
    default:
        std::cout << "Ok it's your dicision\n";
    }
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output << sol.t.size() << "\n";
    output.close();

    plt::named_plot("x(t)", sol.t, sol.x, "-");
    plt::grid(true);
    plt::show();

    return 0;
}