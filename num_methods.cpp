#include "num_methods.h"

matrix numeric_method::euler_method(float_T RelTol, float_T AbsTol)
{
    // Инициализация
    size_t p = 1; // p - порядок метода

    float_T a; // Параметр
    float_T h = 1e-3, h_2; // Шаги
    float_T x, xprev, X; // Координаты
    float_T tmin, tmax, t, tprev; // Время
    float_T eps_cor = 1 / (pow(2, p + 1)), s_coef = 1 / (pow(2, p) - 1), eps, S; // Контроль погрешности

    // Контейнер для результатов решения
    matrix sol;

    // Ввод
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax; // Параметр, начальное значение, интервал интегрирования
    input.close();

    // Расчет
    t = tprev = tmin;
    sol.t.push_back(t); // Сохраняем начальную точку
    sol.x.push_back(x); // Сохраняем начальную точку

    while (t < tmax)
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
            sol.t.push_back(t); // Сохраняем результат
            sol.x.push_back(x); // Сохраняем результат
            if (S < eps * eps_cor) // Погрешность слишком хороша, уменьшаем шаг
                h *= 2;
        }
    }

    if (t > tmax) // Если перешли через границу, выкидываем последнюю точку и считаем прямо на границе
    {
        sol.t.pop_back();
        sol.x.pop_back();
        
        t = tprev;
        x = X = xprev;
        h = tmax - t;

        h_2 = h * 0.5;
        x += h * (a - sin(x));
        X += h_2 * (a - sin(X));
        X += h_2 * (a - sin(X));
        t += h;

        sol.t.push_back(t);
        sol.x.push_back(x);
    }

    return sol;
}

void numeric_method::euler_maruyama_method(matrix& sol, size_t N)
{
    float_T x, a; // Координата, параметр
    float_T sqrt_hd; // Спец. константа
    float_T h, t, tmin, tmax; // Шаг интегрирования, время, интервал интегрирования

    // Создаём генератор случайных чисел со стандартным нормальным распределением
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float_T> norm_dist(0, 1);

    // Ввод данных
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax; // Параметр, начальное значение, интервал интегрирования
    input.close();

    // Расчет
    t = tmin;
    sol.t[0] = t;
    sol.x[0] = x;
    h = (tmax - tmin) / N, sqrt_hd = sqrt(h * D);

    // Вектор случайных величин
    vec_T Z1(N + 1);
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

void numeric_method::hyun_method(matrix& sol, size_t N)
{
    float_T x, a; // Координата, параметр
    float_T sqrt_hd, as; // Спец. константа, спец. выражение
    float_T h, h_2, t, tmin, tmax; // Шаг интегрирования, его половина, время, интервал интегрирования

    // Создаём генератор случайных чисел со стандартным нормальным распределением
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float_T> norm_dist(0, 1);

    // Ввод
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax; // Параметр, начальное значение, интервал интегрирования
    input.close();

    // Расчет
    t = tmin;
    sol.t[0] = t;
    sol.x[0] = x;
    h = (tmax - tmin) / N, h_2 = h * 0.5, sqrt_hd = sqrt(h * D);

    // Вектор случайных величин
    vec_T Z1(N + 1);
    for (size_t i = 1; i <= N; ++i) {
        Z1[i] = sqrt_hd * norm_dist(gen);
    }

    for (size_t i = 1; i <= N; ++i)
    {
        as = a - sin(x);
        x += h_2 * (as + a - sin(x + h * as + Z1[i])) + Z1[i];
        t += h;

        sol.t[i] = t;
        sol.x[i] = x;        
    }
}

void numeric_method::stoch_rk4_method(matrix& sol, size_t N)
{
    float_T x, k1, k2, k3, k4, a; // Координата, параметр
    float_T sqrt_hd, _1_6; // Спец. константа 1, спец. константа 2
    float_T h, h_2, t, tmin, tmax; // Шаг интегрирования, его половина, время, интервал интегрирования

    // Создаём генератор случайных чисел со стандартным нормальным распределением
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float_T> norm_dist(0, 1);

    // Ввод
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax; // Параметр, начальное значение, интервал интегрирования
    input.close();

    // Расчет
    t = tmin;
    sol.t[0] = t;
    sol.x[0] = x;
    h = (tmax - tmin) / N;
    h_2 = h * 0.5, sqrt_hd = sqrt(h * D), _1_6 = 1.0 / 6;

    // Вектор случайных величин
    vec_T Z1(N + 1);
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