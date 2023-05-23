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
    float_T a, x0, tmin, tmax, D; // Интенсивность шума

    // Ввод параметров
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();

    // Ввод параметра, начального значения и диапазона интегрирования
    std::ifstream input("input.txt");
    input >> a >> x0 >> tmin >> tmax; // Параметр, начальное значение, интервал интегрирования
    input.close();
    
    // Начинаем расчёт, засекаем время
    auto start_time = std::chrono::steady_clock::now();

    // Объявляем переменные
    matrix sol;
    sol.t = new float_T[n + 1];
    sol.x = new float_T[n + 1];

    numeric_method solver(D);

    switch (method_number)
    {
    case 1:
        solver.euler_method(sol, n, a, x0, tmin, tmax, 1e-6, 1e-6); // Считаем методом Эйлера, без шума
        break;
    case 2:
        solver.euler_maruyama_method(sol, n, a, x0, tmin, tmax); // Метод Эйлера-Маруямы
        break;
    case 3:
        solver.hyun_method(sol, n, a, x0, tmin, tmax); // Метод Хюна
        break;
    case 4:
        solver.stoch_rk4_method(sol, n, a, x0, tmin, tmax); // Стохастический метод РК4
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
    output.close();


    std::vector<float_T> t(n + 1), x(n + 1);
    for (size_t i = 0; i <= n; ++i)
    {
        t[i] = sol.t[i];
        x[i] = sol.x[i];
    }
    delete[] sol.t;
    delete[] sol.x;

    // Строим график
    plt::named_plot("x(t)", t, x, "-");
    plt::grid(true);
    plt::legend();
    plt::show();

    return 0;
}