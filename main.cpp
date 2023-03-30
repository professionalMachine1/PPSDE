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

    int32_t method_number; // Номер метода
    std::cin >> method_number;

    size_t n; // Количество участков интегрирования
    float_T D; // Интенсивность шума

    // Ввод параметров
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();
    
    // Начинаем расчёт, засекаем время
    auto start_time = std::chrono::steady_clock::now();

    // Объявляем переменные
    matrix sol = { vec_T(n + 1), vec_T(n + 1) };
    numeric_method solver(D);

    switch (method_number)
    {
    case 1:
        sol = matrix();
        sol = solver.euler_method(1e-6, 1e-6); // Считаем методов Эйлера, без шума
        break;
    case 2:
        solver.euler_maruyama_method(sol, n); // Метод Эйлера-Маруямы
        break;
    case 3:
        solver.hyun_method(sol, n); // Метод Хюна
        break;
    case 4:
        solver.stoch_rk4_method(sol, n); // Стохастический метод РК4
        break;
    default:
        std::cout << "Ok it's your dicision\n";
    }

    // Засекаем время завершения расчёта
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output << sol.t.size() << "\n";
    output.close();

    // Строим график
    plt::named_plot("x(t)", sol.t, sol.x, "-");
    plt::grid(true);
    plt::legend();
    plt::show();

    return 0;
}