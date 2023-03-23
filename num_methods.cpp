#include "num_methods.h"

matrix numeric_method::euler_method(float_T eps, float_T h)
{
    // Инициализация
    size_t p = 1; // p - порядок метода

    float_T h_2;
    float_T x, xprev, X, tmin, tmax, t, tprev, a, S;
    float_T eps_cor = eps / (pow(2, p + 1)), s_coef = 1 / (pow(2, p) - 1);

    matrix val_list;

    // Ввод
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax;
    input.close();

    // Расчет
    tprev = t = tmin;
    while (t < tmax)
    {
        tprev = t;
        xprev = x;
        X = x;

        h_2 = h * 0.5;
        x = x + h * (a - sin(x));
        X = X + h_2 * (a - sin(X));
        X = X + h_2 * (a - sin(X));

        S = fabs(X - x) * s_coef;
        if (S > eps)
        {
            x = xprev;
            h *= 0.5;
        }
        else
        {
            t += h;
            val_list.t.push_back(t);
            val_list.x.push_back(x);
            if (S < eps_cor)
                h *= 2;
        }
    }

    if (t > tmax)
    {
        val_list.t.pop_back();
        val_list.x.pop_back();
        
        t = tprev;
        x = X = xprev;
        h = tmax - t;

        h_2 = h * 0.5;
        x = x + h * (a - sin(x));
        X = X + h_2 * (a - sin(X));
        X = X + h_2 * (a - sin(X));
        t += h;

        val_list.t.push_back(t);
        val_list.x.push_back(x);
    }

    return val_list;
}

matrix numeric_method::stoch_euler_method(size_t N)
{
    float_T Z1, h, x, tmin, tmax, t, a;

    matrix val_list = { std::vector<float_T>(N + 1), 
        std::vector<float_T>(N + 1) };

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float_T> norm_dist(0, 1);

    // Ввод данных
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax;
    input.close();

    // Расчет
    t = tmin;
    val_list.t[0] = t;
    val_list.x[0] = x;

    h = (tmax - tmin) / N;
    for (size_t i = 1; i <= N; ++i)
    {
        Z1 = sqrt(h * D) * norm_dist(gen);
        x = x + h * (a - sin(x)) + Z1;
        t += h;

        val_list.t[i] = t;
        val_list.x[i] = x;
    }

    return val_list;
}

matrix numeric_method::stoch_hyun_method(size_t N)
{
    float_T Z1, x, tmin, tmax, t, a, h, h_2;

    matrix val_list = { std::vector<float_T>(N + 1), 
        std::vector<float_T>(N + 1) };

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float_T> norm_dist(0, 1);

    // Ввод
    std::ifstream input("input.txt");
    input >> a >> x >> tmin >> tmax;
    input.close();

    // Расчет
    t = tmin;
    val_list.t[0] = t;
    val_list.x[0] = x;
    h = (tmax - tmin) / N, h_2 = h * 0.5;

    for (size_t i = 0; i <= N; ++i)
    {
        Z1 = sqrt(h * D) * norm_dist(gen);
        x = x + h_2 * (a - sin(x) + a - sin(x + h * (a - sin(x)) + Z1)) + Z1;
        t += h;

        val_list.t[i] = t;
        val_list.x[i] = x;        
    }

    return val_list;
}
