# PPSDE
Решается задача интегрирования стохастического дифференциального уравнения

\frac{dx}{dt} = a - \sin(x) + \xi(t), \xi ~ \N(0, D)

Задача решается несколькими способами:
  * Методом Эйлера с контролем локальной погрешности, случайная величина \xi не учитывается
  * Стохастическим методо Эйлера
  * Методом Хюна