import numpy as np
import matplotlib.pyplot as plt


def rof(module, alpha, c_coef):
    return c_coef * module / (1 - np.sin(alpha))

def z_min(ha_coef, alpha):                           # Расчет минимального количества зубьев
    return 2 * ha_coef / np.power(np.sin(alpha), 2)

def rad_del(module, z):                              # Расчет делительного радиуса
    return module * z / 2

def rad_osn(module, z, alpha):                       # Расчет основного радиуса
    return rad_del(module, z) * np.cos(alpha)

def x_min(z, zmin, ha_coef):                         # Расчет минимального смещения
    return ha_coef * (zmin - z) / zmin

def angle_inv(a):                                    # Поиск угла инволюты по методу Ньютона
    x = np.pi / 4
    while np.fabs(np.tan(x) - x - a) > 1e-4:
        x = x - (np.tan(x) - x - a) / pow((np.tan(x)), 2)
    return x

def alpha_w(x1, x2, z1, z2, alpha):                  # Угол зацепления
    return angle_inv(np.tan(alpha) - alpha + 2 * np.tan(alpha) * (x1 + x2) / (z1 + z2))

def y(x1, x2, z1, z2, alpha):                        # Коэффицент воспринимаемого смещения
    return (z1 + z2) * (np.cos(alpha) / np.cos(alpha_w(x1, x2, z1, z2, alpha)) - 1) / 2

def dy(x1, x2, y):                                   # Коэффицент уравнительного смещения
    return x1 + x2 - y

def interaxis(r1, r2):                               # Межосевое расстояние
    return r1 + r2

def ra(module, z, x, dy, ha_coef):                   # Радиус окружности вершин
    return module * (z / 2 + ha_coef + x - dy)

def rf(del_rad, x, module, c_coef, ha_coef):         # Радиус окружности впадин
    return del_rad + x * module - ha_coef * module - c_coef * module

def height(module, dy, ha_coef, c_coef):             # Высота зуба
    return module * (2 * ha_coef + c_coef - dy)

def S(module, x, alpha):                             # Толщина зуба по делительной окружности
    return module * (np.pi / 2 + 2 * x * np.tan(alpha))

def Sa(ra, S, r, rb, alpha):                         # Толщина зуба по окружности вершин
    return 2 * ra * (S / (2 * r) + np.tan(alpha) - alpha - np.tan(np.arccos(rb / ra)) + np.arccos(rb / ra))

def eps_alpha(z1, z2, ra1, ra2, rb1, rb2, alpha_w):  # Коэффицент перекрытия прямозубой передачи
    return z1 / (2 * np.pi) * (np.tan(np.arccos(rb1 / ra1)) - np.tan(alpha_w)) + z2 / (2 * np.pi) * (np.tan(np.arccos(rb2 / ra2)) - np.tan(alpha_w))

def lambda_f(z1, z2, ra, rb, alpha_w):               # Коэффицент удельного скольжения
    if z1 > z2:
        return z2 * (1 + z2 / z1) * (np.tan(np.arccos(rb / ra)) - np.tan(alpha_w)) / ((z1 + z2) * np.tan(alpha_w) - z2 * np.tan(np.arccos(rb / ra)))
    else:
        return z2 * (1 + z1 / z2) * (np.tan(np.arccos(rb / ra)) - np.tan(alpha_w)) / ((z1 + z2) * np.tan(alpha_w) - z2 * np.tan(np.arccos(rb / ra)))

def thetta(z1, z2, alpha_w, alpha):                  # Определение коэффицента удельного давления
    return 2 * (z1 + z2) / (z1 * z2 * np.tan(alpha_w) * np.cos(alpha))

def graph(start, end, step, **kwargs):               # Получение массива точек графиков

    zmin = z_min(kwargs["ha_coef"], kwargs["alpha"])
    border = {"START": start, "END": end, "STEP": step}
    x_min_1 = x_min(kwargs["z1"], zmin, kwargs["ha_coef"])
    x_min_2 = x_min(kwargs["z2"], zmin, kwargs["ha_coef"])
    interval = int((border["END"] - border["START"]) / border["STEP"])

    #Графики eps_alpha, [Sa1], [Sa2], lambda1 = lambda2

    x1 = np.arange(x_min_1, border["END"], border["STEP"])
    x2 = np.arange(x_min_2, border["END"], border["STEP"])

    x_res_sa1 = np.array([])
    y_res_sa1 = np.array([])
    x_res_sa2 = np.array([])
    y_res_sa2 = np.array([])
    x_res_epsa = np.array([])
    y_res_epsa = np.array([])
    x_res_lambda = np.array([])
    y_res_lambda = np.array([])
    x_res_thetta = np.array([])
    y_res_thetta = np.array([])

    for i in x1:
        for j in x2:
            
            y_t = y(i, j, kwargs["z1"], kwargs["z2"], kwargs["alpha"])
            dy_t = dy(i, j, y_t)
            alpha_w_t = alpha_w(i, j, kwargs["z1"], kwargs["z2"], kwargs["alpha"])
            ra1 = ra(kwargs["module"], kwargs["z1"], i, dy_t, kwargs["ha_coef"])
            ra2 = ra(kwargs["module"], kwargs["z2"], j, dy_t, kwargs["ha_coef"])
            S1 = S(kwargs["module"], i, kwargs["alpha"])
            S2 = S(kwargs["module"], j, kwargs["alpha"])
            r1 = rad_del(kwargs["module"], kwargs["z1"])
            r2 = rad_del(kwargs["module"], kwargs["z2"])
            rb1 = rad_osn(kwargs["module"], kwargs["z1"], kwargs["alpha"])
            rb2 = rad_osn(kwargs["module"], kwargs["z2"], kwargs["alpha"])
            Sa1 = Sa(ra1, S1, r1, rb1, kwargs["alpha"])
            Sa2 = Sa(ra2, S2, r2, rb2, kwargs["alpha"])
            eps_alpha_t = eps_alpha(kwargs["z1"], kwargs["z2"], ra1, ra2, rb1, rb2, alpha_w_t)
            lambda1 = lambda_f(z1=kwargs["z1"], z2=kwargs["z2"], ra=ra2, rb=rb2, alpha_w=alpha_w_t)
            lambda2 = lambda_f(z1=kwargs["z2"], z2=kwargs["z1"], ra=ra1, rb=rb1, alpha_w=alpha_w_t)
            thetta_t = thetta(kwargs["z1"], kwargs["z2"], alpha_w_t, kwargs["alpha"])

            if np.fabs((Sa1 - kwargs["Sa_coef"] * kwargs["module"])) < kwargs["epsilon_sa"]:
                x_res_sa1 = np.append(x_res_sa1, i)
                y_res_sa1 = np.append(y_res_sa1, j)
            if np.fabs((Sa2 - kwargs["Sa_coef"] * kwargs["module"])) < kwargs["epsilon_sa"]:
                x_res_sa2 = np.append(x_res_sa2, i)
                y_res_sa2 = np.append(y_res_sa2, j)
            if np.fabs((eps_alpha_t - kwargs["eps_alpha_coef"])) < kwargs["epsilon_epsa"]:
                x_res_epsa = np.append(x_res_epsa, i)
                y_res_epsa = np.append(y_res_epsa, j)
            if np.fabs(lambda1 - lambda2) < kwargs["epsilon_lambda"] and lambda1 > 0 and lambda2 > 0:
                x_res_lambda = np.append(x_res_lambda, i)
                y_res_lambda = np.append(y_res_lambda, j)
            if np.fabs(thetta_t - kwargs["thetta_max"]) < kwargs["epsilon_thetta"]:
                x_res_thetta = np.append(x_res_thetta, i)
                y_res_thetta = np.append(y_res_thetta, j)


    #График x_min_1

    xmin1 = np.full(interval, x_min_1)

    #График x_min_2

    xmin2 = np.full(interval, x_min_2)

    #Зарисовка графиков

    plt.title("Блокирующий контур")
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.grid()
    plt.plot(xmin1, x2, "red")
    plt.plot(x1, xmin2, "red")
    plt.plot(x_res_sa1, y_res_sa1, "green")
    plt.plot(x_res_sa2, y_res_sa2, "green")
    plt.plot(x_res_epsa, y_res_epsa, "purple")
    plt.plot(x_res_lambda, y_res_lambda, "brown")
    plt.plot(x_res_thetta, y_res_thetta, "yellow")

    print("Значения с прямой lambda1 = lambda2:")
    for i, j in zip(x_res_lambda, y_res_lambda):
        print(f"x1: {i:.4}, x2: {j:.4}")    

def analyze_single_condition(q1, q2, **kwargs):      # Определение параметров при фиксированном x1, x2
    
    rof_t = rof(kwargs["module"], kwargs["alpha"], kwargs["c_coef"])
    zmin = z_min(kwargs["ha_coef"], kwargs["alpha"])
    minx1 = x_min(kwargs["z1"], zmin, kwargs["ha_coef"])
    minx2 = x_min(kwargs["z2"], zmin, kwargs["ha_coef"])
    y_t = y(q1, q2, kwargs["z1"], kwargs["z2"], kwargs["alpha"])
    dy_t = dy(q1, q2, y_t)
    h = height(kwargs["module"], dy_t, kwargs["ha_coef"], kwargs["c_coef"])
    alpha_w_t = alpha_w(q1, q2, kwargs["z1"], kwargs["z2"], kwargs["alpha"])
    S1 = S(kwargs["module"], q1, kwargs["alpha"])
    S2 = S(kwargs["module"], q2, kwargs["alpha"])
    r1 = rad_del(kwargs["module"], kwargs["z1"])
    r2 = rad_del(kwargs["module"], kwargs["z2"])
    rf1 = rf(r1, q1, kwargs["module"], kwargs["c_coef"], kwargs["ha_coef"])
    rf2 = rf(r2, q2, kwargs["module"], kwargs["c_coef"], kwargs["ha_coef"])
    rb1 = rad_osn(kwargs["module"], kwargs["z1"], kwargs["alpha"])
    rb2 = rad_osn(kwargs["module"], kwargs["z2"], kwargs["alpha"])
    ra1 = ra(kwargs["module"], kwargs["z1"], q1, dy_t, kwargs["ha_coef"])
    ra2 = ra(kwargs["module"], kwargs["z2"], q2, dy_t, kwargs["ha_coef"])
    Sa1 = Sa(ra1, S1, r1, rb1, kwargs["alpha"])
    Sa2 = Sa(ra2, S2, r2, rb2, kwargs["alpha"])
    eps_alpha_t = eps_alpha(kwargs["z1"], kwargs["z2"], ra1, ra2, rb1, rb2, alpha_w_t)
    lambda1 = lambda_f(z1=kwargs["z1"], z2=kwargs["z2"], ra=ra2, rb=rb2, alpha_w=alpha_w_t)
    lambda2 = lambda_f(z1=kwargs["z2"], z2=kwargs["z1"], ra=ra1, rb=rb1, alpha_w=alpha_w_t)
    thetta_t = thetta(kwargs["z1"], kwargs["z2"], alpha_w_t, kwargs["alpha"])

    print(f"Минимальное количество зубьев:                                    {int(zmin):> 5}",
          f"Минимальное смещение для первого колеса с отсутствием подрезания: {minx1:> 5.3f} мм",
          f"Минимальное смещение для второго колеса с отсутствием подрезания: {minx2:> 5.3f} мм",
          f"Угол зацепления:                                                  {alpha_w_t * 180/ np.pi:> 5.3f} DEG",
          f"Уравнительное смещение:                                           {dy_t:> 5.3f} мм",
          f"Воспринимаемое смещение:                                          {y_t:> 5.3f} мм",
          f"Радиус первого колеса по окружности вершин:                       {ra1:> 5.3f} мм",
          f"Радиус второго колеса по окружности вершин:                       {ra2:> 5.3f} мм",
          f"Радиус первого колеса по окружности впадин:                       {rf1:> 5.3f} мм",
          f"Радиус первого колеса по окружности впадин:                       {rf2:> 5.3f} мм",
          f"Высота зуба:                                                      {h:> 5.3f} мм",
          f"Толщина зуба первого колеса по делительной окружности:            {S1:> 5.3f} мм",
          f"Толщина зуба второго колеса по делительной окружности:            {S2:> 5.3f} мм",
          f"Толщина зуба первого колеса по окружности вершин:                 {Sa1:> 5.3f} мм",
          f"Толщина зуба второго колеса по окружности вершин:                 {Sa2:> 5.3f} мм",
          f"Коэффицент торцевого перекрытия:                                  {eps_alpha_t:> 5.3f}",
          f"Коэффицент скольжения для первого колеса:                         {lambda1:> 5.3f}",
          f"Коэффицент скольжения для второго колеса:                         {lambda2:> 5.3f}",
          f"Коэффицент удельного давления:                                    {thetta_t:> 5.3f} DEG", sep="\n")

def planet(**kwargs):                                # Определение параметров планетарного механизма

    res = []
    for z20 in range(kwargs["zmin_extern"], kwargs["zmax"]):
        for z21 in range(kwargs["zmin_extern"], kwargs["zmax"]):
            for z22 in range(kwargs["zmin_intern"], kwargs["zmax"]):
                if (z22 - z21 >= 8 and np.sqrt(3)/2 > (z21 + 2 * kwargs["ha_coef"]) / (z20 + z21) and
                    z20 + z21 == z22 - z21 and 1 + z22 / z20 == kwargs["u20h"]):
                    res.append((z20, z21, z22))

    print(*res, sep="\n")

    ans = [int(i) for i in input("Введите через пробел количество зубьев у 21, 22 и 23 колес: ").split()]
    res = []

    for i in range(100):
        for j in range(100):
            if ans[0] * kwargs["u20h"] * (1 + 3 * i) == j:
                res.append((i, j))
                break

    print(f"При П == {res[0][0]} и Z == {res[0][1]} сборка возможна.\n"
        f"Угол поворота водила fi_h == {(2 * np.pi / 3 + 2 * np.pi * res[0][0]) * 180 / np.pi}")