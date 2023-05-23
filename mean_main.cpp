#include "num_methods.h"
#include "matplotlibcpp.h"
#include <chrono>
#include <omp.h>

namespace plt = matplotlibcpp;

void vec_sum(vec_T change, const vec_T add, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        change[i] += add[i];
    }
}

void vec_mul_num(vec_T change, const float_T val, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        change[i] *= val;
    }
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t n; // Количество участков интегрирования
    float_T a, x0, tmin, tmax, D, h; // Интенсивность шума

    // Ввод параметров
    std::ifstream params("params.txt");
    params >> n >> D;
    params.close();

    // Ввод параметра, начального значения и диапазона интегрирования
    std::ifstream input("input.txt");
    input >> a >> x0 >> tmin >> tmax; // Параметр, начальное значение, интервал интегрирования
    input.close();

    size_t N; // Количество реализаций
    std::cout << "Enter number of implementations\n";
    std::cin >> N;

    // Начинаем расчёт и засекаем время
    auto start_time = std::chrono::steady_clock::now();

    // Объявляем переменные
    matrix cum_sol;
    cum_sol.t = new float_T[n + 1];
    cum_sol.x = new float_T[n + 1];

    cum_sol.t[0] = tmin;
    h = (tmax - tmin) / (1.0 * n);
    for (size_t i = 1; i <= n; ++i) {
        cum_sol.t[i] = cum_sol.t[i - 1] + h;
    }

    // Получаем реализации и суммируем их
    // Последовательная версия
    // numeric_method solver(D);

    // matrix sol;
    // sol.t = new float_T[n + 1];
    // sol.x = new float_T[n + 1];

    // for (size_t i = 0; i < N; ++i) 
    // {
    //     solver.hyun_method(sol, n, a, x0, tmin, tmax);
    //     vec_sum(cum_sol.x, sol.x, n + 1);
    // }
    // // Строим решение без шума для проверки сходимости
    // solver.euler_method(sol, n, a, x0, tmin, tmax, 1e-6, 1e-6);
    
    // Параллельная версия
    #pragma omp parallel
    {
        vec_T cur_cum_sol = new float_T[n + 1];
        numeric_method solver(D);
        matrix sol;

        sol.t = new float_T[n + 1];
        sol.x = new float_T[n + 1];

        #pragma omp for
        for (int32_t i = 0; i < N; ++i) 
        {
            solver.hyun_method(sol, n, a, x0, tmin, tmax);
            vec_sum(cur_cum_sol, sol.x, n + 1);
        }
        #pragma omp critical
        {
            vec_sum(cum_sol.x, cur_cum_sol, n + 1);
        }
        delete[] sol.t;
        delete[] sol.x;
        delete[] cur_cum_sol;
    }
    // Строим решение без шума для проверки сходимости
    matrix sol;
    sol.t = new float_T[n + 1];
    sol.x = new float_T[n + 1];

    numeric_method solver(D);
    solver.euler_method(sol, n, a, x0, tmin, tmax, 1e-6, 1e-6);

    // Усредняем реализации
    vec_mul_num(cum_sol.x, (1.0 / N), n + 1);

    // Засекаем время окончания расчёта
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output.close();

    std::vector<float_T> tv(n + 1), xv(n + 1), tr(n + 1), xr(n + 1);
    for (size_t i = 0; i <= n; ++i)
    {
        tv[i] = cum_sol.t[i];
        xv[i] = cum_sol.x[i];
        tr[i] = sol.t[i];
        xr[i] = sol.x[i];
    }
    delete[] cum_sol.t;
    delete[] cum_sol.x;
    delete[] sol.t;
    delete[] sol.x;

    // Строим график усреднённых реализаций
    plt::named_plot("~x(t)", tv, xv, "-");
    plt::named_plot("x(t)", tr, xr, "-");
    plt::grid(true);
    plt::legend();
    plt::show();

    return 0;
}