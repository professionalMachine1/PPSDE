#include "num_methods.h"
#include "matplotlibcpp.h"
#include <chrono>

namespace plt = matplotlibcpp;

vec_T& operator+=(vec_T& change, const vec_T& add)
{
    size_t n = change.size();
    for (size_t i = 0; i < n; ++i) {
        change[i] += add[i];
    }

    return change;
}

vec_T& operator*=(vec_T& change, const float_T val)
{
    size_t n = change.size();
    for (size_t i = 0; i < n; ++i) {
        change[i] *= val;
    }
    return change;
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t N; // Количество реализаций
    std::cout << "Enter number of implementations\n";
    std::cin >> N;

    size_t n; // Количество участков интегрирования
    float_T D; // Интенсивность шума

    // Ввод параметров
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();

    // Левая и правая границы ямки
    float_T lb = -2 * asin(1), rb = 2 * asin(1);

    // Начинаем расчёт, засекаем время
    auto start_time = std::chrono::steady_clock::now();

    // Объявляем переменные
    matrix sol = { vec_T(n + 1), vec_T(n + 1) }, density = sol;
    numeric_method solver(D);
    
    // Получаем реализации и суммируем их
    // Последовательная версия
    for (size_t i = 0; i < N; ++i)
    {
        solver.hyun_method(sol, n);
        for (size_t j = 0; j <= n; ++j) { // Считаем попадания точек в интервал
            if (sol.x[j] > lb && sol.x[j] < rb)
                density.x[j] += 1;
        }
    }
    
    // Параллельная версия
    // #pragma omp parallel private(sol, solver)
    // {
    //     vec_T cur_density(n + 1);
    //     #pragma omp for
    //     for (size_t i = 0; i < N; ++i)
    //     {
    //         solver.hyun_method(sol, n);
    //         for (size_t j = 0; j <= n; ++j) { // Считаем попадания точек в интервал
    //             if (sol.x[j] > lb && sol.x[j] < rb)
    //                 cur_density[j] += 1;
    //         }
    //     }
    //     #pragma omp critical
    //     density.x += cur_density;
    // }

    density.t = sol.t; // Сохраняем время для распределения
    density.x *= (1.0 / N); // Усредняем вхождения

    // Засекаем время окончания расчёта
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output << sol.t.size() << "\n";
    output.close();

    // Строим график
    plt::named_plot("p(t)", density.t, density.x, "-");
    plt::grid(true);
    plt::legend();
    plt::show();

    return 0;
}