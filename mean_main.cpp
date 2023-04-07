#include "num_methods.h"
#include "matplotlibcpp.h"
#include <chrono>
#include <omp.h>

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

    size_t n; // Количество участков интегрирования
    float_T D; // Интенсивность шума

    // Ввод параметров
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();

    size_t N; // Количество реализаций
    std::cout << "Enter number of implementations\n";
    std::cin >> N;
    
    // Начинаем расчёт и засекаем время
    auto start_time = std::chrono::steady_clock::now();

    // Объявляем переменные
    matrix cum_sol = { vec_T(n + 1), vec_T(n + 1) };

    // Получаем реализации и суммируем их
    // Последовательная версия
    // matrix sol = { vec_T(n + 1), vec_T(n + 1) };
    // numeric_method solver(D);

    // for (size_t i = 0; i < N; ++i) 
    // {
    //     solver.hyun_method(sol, n);
    //     cum_sol.x += sol.x;
    // }
    // cum_sol.t = sol.t; // Устанавливаем время в общее решение
    // // Строим решение без шума для проверки сходимости
    // sol = matrix();
    // sol = solver.euler_method(1e-6, 1e-6);
    
    // Параллельная версия
    #pragma omp parallel
    {
        vec_T cur_cum_sol(n + 1);
        numeric_method solver(D);
        matrix sol = { vec_T(n + 1), vec_T(n + 1) };
        #pragma omp for
        for (int32_t i = 0; i < N; ++i) 
        {
            solver.hyun_method(sol, n);
            cur_cum_sol += sol.x;
        }
        #pragma omp critical
        cum_sol.x += cur_cum_sol;
        #pragma omp single
        cum_sol.t = sol.t;
    }
    // Строим решение без шума для проверки сходимости
    matrix sol = matrix();
    numeric_method solver(D);
    sol = solver.euler_method(1e-6, 1e-6);

    cum_sol.x *= (1.0 / N); // Усредняем реализации

    // Засекаем время окончания расчёиа
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output << cum_sol.t.size() << "\n";
    output.close();

    // Строим график усреднённых реализаций
    plt::named_plot("~x(t)", cum_sol.t, cum_sol.x, "-");
    plt::named_plot("x(t)", sol.t, sol.x, "-");
    plt::grid(true);
    plt::legend();
    plt::show();

    return 0;
}