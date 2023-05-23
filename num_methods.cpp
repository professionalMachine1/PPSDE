#include "num_methods.h"

void numeric_method::euler_method(matrix& sol, size_t N, float_T a, float_T x0, 
        float_T tmin, float_T tmax, float_T RelTol, float_T AbsTol)
{
    // Инициализация
    size_t p = 1; // p - порядок метода

    float_T t, tprev; // Время
    float_T x, xprev, X; // Координаты
    float_T h = (tmax - tmin) / (1.0 * N), h_2; // Шаги
    float_T eps_cor = 1 / (pow(2, p + 1)), s_coef = 1 / (pow(2, p) - 1), eps, S; // Контроль погрешности

    sol.x[0] = x0;
    sol.t[0] = tmin;
    for (size_t i = 1; i <= N; ++i) {
        sol.t[i] = sol.t[i - 1] + h;
    }

    // Расчет
    x = x0;
    t = tmin;
    for (size_t i = 1; i <= N; )
    {
        // Запоминаем текущие значения
        X = x;
        xprev = x;
        tprev = t;

        // Просчитываем следующую точку
        h_2 = h * 0.5;
        x += h * (a - sin(x));
        X += h_2 * (a - sin(X));
        X += h_2 * (a - sin(X));

        // Определяем погрешность на шаге
        S = fabs(X - x) * s_coef;
        eps = std::max(RelTol * abs(x), AbsTol);

        if (S > eps) // Погрешность велика, возвращаемся к предыдущей точке, уменьшаем шаг
        {
            x = xprev;
            h *= 0.5;
        }
        else // Погрешность хороша, сохраняем решение
        {
            t += h; // Увеличиваем время
            if (t >= sol.t[i]) { // Сохраняем результат
                sol.x[i] = xprev + (sol.t[i] - tprev) * (a - sin(xprev));
                // std::cout << sol.x[i] << "\n";
                ++i;
            }
            if (S < eps * eps_cor) // Погрешность слишком хороша, уменьшаем шаг
                h *= 2;
        }
    }
}

void numeric_method::euler_maruyama_method(matrix& sol, size_t N, float_T a, 
    float_T x0, float_T tmin, float_T tmax)
{
    float_T x = x0; // Координата, параметр
    float_T sqrt_hd; // Спец. константа
    float_T h, t; // Шаг интегрирования, время, интервал интегрирования

    // Создаём генератор случайных чисел со стандартным нормальным распределением
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    boost::normal_distribution<float_T> norm_dist(0, 1);

    // Расчет
    t = tmin;
    sol.t[0] = t;
    sol.x[0] = x;
    h = (tmax - tmin) / N, sqrt_hd = sqrt(h * D);

    // Вектор случайных величин
    vec_T Z1 = new float_T[N + 1];
    for (size_t i = 1; i <= N; ++i) {
        Z1[i] = sqrt_hd * norm_dist(gen);
    }

    for (size_t i = 1; i <= N; ++i)
    {
        x += h * (a - sin(x)) + Z1[i];
        t += h;

        sol.t[i] = t;
        sol.x[i] = x;
    }
}

void numeric_method::hyun_method(matrix& sol, size_t N, float_T a, 
    float_T x0, float_T tmin, float_T tmax)
{
    float_T x = x0; // Координата, параметр
    float_T sqrt_hd, as; // Спец. константа, спец. выражение
    float_T h, h_2, t; // Шаг интегрирования, его половина, время, интервал интегрирования

    // Создаём генератор случайных чисел со стандартным нормальным распределением
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    boost::normal_distribution<float_T> norm_dist(0, 1);

    // Расчет
    t = tmin;
    sol.t[0] = t;
    sol.x[0] = x;
    h = (tmax - tmin) / N, h_2 = h * 0.5, sqrt_hd = sqrt(h * D);

    // Вектор случайных величин
    vec_T Z1 = new float_T[N + 1];
    for (size_t i = 1; i <= N; ++i) {
        Z1[i] = sqrt_hd * norm_dist(gen);
    }

    for (size_t i = 1; i <= N; ++i)
    {
        as = a - sin(x);
        // x += h_2 * (as + a - sin(x + Z1[i])) / (1 - h_2 * cos(x)) + Z1[i];
        x += h_2 * (as + a - sin(x + h * as + Z1[i])) + Z1[i];
        t += h;

        sol.t[i] = t;
        sol.x[i] = x;        
    }

    delete[] Z1;
}

void numeric_method::stoch_rk4_method(matrix& sol, size_t N, float_T a, 
    float_T x0, float_T tmin, float_T tmax)
{
    float_T x = x0, k1, k2, k3, k4; // Координата, параметр
    float_T sqrt_hd, _1_6; // Спец. константа 1, спец. константа 2
    float_T h, h_2, t; // Шаг интегрирования, его половина, время, интервал интегрирования

    // Создаём генератор случайных чисел со стандартным нормальным распределением
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    boost::normal_distribution<float_T> norm_dist(0, 1);

    // Расчет
    t = tmin;
    sol.t[0] = t;
    sol.x[0] = x;
    h = (tmax - tmin) / N;
    h_2 = h * 0.5, sqrt_hd = sqrt(h * D), _1_6 = 1.0 / 6;

    // Вектор случайных величин
    vec_T Z1 = new float_T[N + 1];
    for (size_t i = 1; i <= N; ++i) {
        Z1[i] = sqrt_hd * norm_dist(gen);
    }

    for (size_t i = 1; i <= N; ++i)
    {
        k1 = (a - sin(x)) * h + Z1[i];
        k2 = (a - sin(x + k1 * 0.5)) * h + Z1[i];
        k3 = (a - sin(x + k2 * 0.5)) * h + Z1[i];
        k4 = (a - sin(x + k3)) * h + Z1[i];
        x += (k1 + 2 * k2 + 2 * k3 + k4) * _1_6;
        t += h;

        sol.t[i] = t;
        sol.x[i] = x;        
    }
}
