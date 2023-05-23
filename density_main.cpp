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

    size_t N; // Количество реализаций
    std::cout << "Enter number of implementations\n";
    std::cin >> N;

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

    // Левая и правая границы ямки
    float_T lb = -2 * asin(1), rb = 2 * asin(1);

    // Начинаем расчёт, засекаем время
    auto start_time = std::chrono::steady_clock::now();

    // Объявляем переменную
    matrix density;
    density.t = new float_T[n + 1];
    density.x = new float_T[n + 1];

    density.x[0] = 0.0;
    density.t[0] = tmin;
    h = (tmax - tmin) / (1.0 * n);
    for (size_t i = 1; i <= n; ++i) 
    {
        density.x[i] = 0.0;
        density.t[i] = density.t[i - 1] + h;
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
    //     for (size_t j = 0; j <= n; ++j) { // Считаем попадания точек в интервал
    //         if (sol.x[j] > lb && sol.x[j] < rb)
    //             density.x[j] += 1;
    //     }
    // }

    // Параллельная версия
    #pragma omp parallel
    {
        vec_T cur_density = new float_T[n + 1];
        numeric_method solver(D);

        matrix sol;
        sol.t = new float_T[n + 1];
        sol.x = new float_T[n + 1];

        for (size_t i = 0; i <= n; ++i) {
            cur_density[i] = 0.0;
        }

        #pragma omp for
        for (int32_t i = 0; i < N; ++i)
        {
            solver.hyun_method(sol, n, a, x0, tmin, tmax);
            for (size_t j = 0; j <= n; ++j) { // Считаем попадания точек в интервал
                if (lb < sol.x[j] && sol.x[j] < rb)
                    cur_density[j] += 1;
            }
        }
        #pragma omp critical
        {
            vec_sum(density.x, cur_density, n + 1);
        }
        delete[] sol.t;
        delete[] sol.x;
        delete[] cur_density;
    }

    // Усредняем вхождения
    vec_mul_num(density.x, (1.0 / N), n + 1);

    // Засекаем время окончания расчёта
    auto end_time = std::chrono::steady_clock::now();

    // Вывод
    std::ofstream output("output.txt");
    output << "solve time = " << std::chrono::duration_cast<
        std::chrono::milliseconds>(end_time - start_time).count() << "\n";
    output.close();

    std::vector<float_T> t(n + 1), x(n + 1);
    for (size_t i = 0; i <= n; ++i)
    {
        t[i] = density.t[i];
        x[i] = density.x[i];
    }
    delete[] density.t;
    delete[] density.x;

    // Строим график
    plt::named_semilogy("p(t)", t, x, "-");
    plt::grid(true);
    plt::legend();
    plt::show();

    return 0;
}