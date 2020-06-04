# Cubic splines
import time

from scipy import misc

from polynomes import *

import polynomes


# linear complexity
def cubic_splines(f, a, b, n):
    xcoords = [a + k * (b - a) / n for k in range(n + 1)]
    ycoords = [f(k) for k in xcoords]
    dycoords = [misc.derivative(f, k) for k in xcoords]
    for k in range(n):
        a0 = ycoords[k + 1] / (xcoords[k + 1] - xcoords[k])
        a1 = ycoords[k] / (xcoords[k] - xcoords[k + 1])
        a2 = (dycoords[k] - a0 - a1) / ((xcoords[k] - xcoords[k + 1]) ** 2)
        a3 = (dycoords[k + 1] - a0 - a1) / ((xcoords[k + 1] - xcoords[k]) ** 2)
        nu0 = -a0 * xcoords[k] - a1 * xcoords[k + 1] - a2 * xcoords[k] * xcoords[k + 1] ** 2 - a3 * xcoords[k] ** 2 * \
              xcoords[k + 1]
        nu1 = a0 + a1 + a2 * (xcoords[k + 1] ** 2 + 2 * xcoords[k + 1] * xcoords[k]) + a3 * (
                xcoords[k] ** 2 + 2 * xcoords[k] * xcoords[k + 1])
        nu2 = -a2 * (xcoords[k] + 2 * xcoords[k + 1]) - a3 * (xcoords[k + 1] + 2 * xcoords[k])
        nu3 = a2 + a3
        P = [nu0, nu1, nu2, nu3]
        polynomial_graph(P, xcoords[k], xcoords[k + 1])
    coordsx = np.linspace(a, b, 1000)
    coordsy = [f(x) for x in coordsx]
    plt.plot(coordsx, coordsy, "b-", label="original function")
    plt.legend(loc="best")


# same complexity
def cubic_splines2(f, a, b, n):
    xcoords = [a + k * (b - a) / n for k in range(n + 1)]
    ycoords = [f(k) for k in xcoords]
    dycoords = [misc.derivative(f, k) for k in xcoords]
    for k in range(n):
        M = np.array([(1, xcoords[k], xcoords[k] ** 2, xcoords[k] ** 3),
                      (1, xcoords[k + 1], xcoords[k + 1] ** 2, xcoords[k + 1] ** 3),
                      (0, 1, 2 * xcoords[k], 3 * xcoords[k] ** 2), (0, 1, 2 * xcoords[k + 1], 3 * xcoords[k + 1] ** 2)])
        N = np.array([(ycoords[k]), (ycoords[k + 1]), (dycoords[k]), (dycoords[k + 1])])
        O = np.dot(np.linalg.inv(M), N)
        a0, a1, a2, a3 = O[0], O[1], O[2], O[3]
        P = [a0, a1, a2, a3]
        polynomial_graph(P, xcoords[k], xcoords[k + 1])
    coordsx = np.linspace(a, b, 1000)
    coordsy = [f(x) for x in coordsx]
    plt.plot(coordsx, coordsy, "b-", label="original function")
    plt.legend(loc="best")


# Integral approximation with cubic splines
# linear complexity
def interp_splines(f, a, b, n):
    I = 0
    xcoords = [a + k * (b - a) / n for k in range(n + 1)]
    ycoords = [f(k) for k in xcoords]
    dycoords = [misc.derivative(f, k) for k in xcoords]
    for k in range(n):
        a0 = ycoords[k + 1] / (xcoords[k + 1] - xcoords[k])
        a1 = ycoords[k] / (xcoords[k] - xcoords[k + 1])
        a2 = (dycoords[k] - a0 - a1) / ((xcoords[k] - xcoords[k + 1]) ** 2)
        a3 = (dycoords[k + 1] - a0 - a1) / ((xcoords[k + 1] - xcoords[k]) ** 2)
        nu0 = -a0 * xcoords[k] - a1 * xcoords[k + 1] - a2 * xcoords[k] * xcoords[k + 1] ** 2 - a3 * xcoords[k] ** 2 * \
              xcoords[k + 1]
        nu1 = a0 + a1 + a2 * (xcoords[k + 1] ** 2 + 2 * xcoords[k + 1] * xcoords[k]) + a3 * (
                xcoords[k] ** 2 + 2 * xcoords[k] * xcoords[k + 1])
        nu2 = -a2 * (xcoords[k] + 2 * xcoords[k + 1]) - a3 * (xcoords[k + 1] + 2 * xcoords[k])
        nu3 = a2 + a3
        P = [nu0, nu1, nu2, nu3]
        I += polynomial_integral(P, xcoords[k], xcoords[k + 1])
    return I


# Comparing the temporal efficiency of both methods on common functions
def cubic_splines_test(a, b, n):
    res_f_ref = [np.sin, np.cos, np.tan, np.exp, np.cosh, np.sinh, np.tanh]
    time1, time2 = [], []
    p = len(res_f_ref)
    for i in range(p):
        # method 1
        f = res_f_ref[i]
        start1 = time.time()
        xcoords = [a + k * (b - a) / n for k in range(n + 1)]
        ycoords = [f(k) for k in xcoords]
        dycoords = [misc.derivative(f, k) for k in xcoords]
        for k in range(n):
            a0 = ycoords[k + 1] / (xcoords[k + 1] - xcoords[k])
            a1 = ycoords[k] / (xcoords[k] - xcoords[k + 1])
            a2 = (dycoords[k] - a0 - a1) / ((xcoords[k] - xcoords[k + 1]) ** 2)
            a3 = (dycoords[k + 1] - a0 - a1) / ((xcoords[k + 1] - xcoords[k]) ** 2)
            nu0 = -a0 * xcoords[k] - a1 * xcoords[k + 1] - a2 * xcoords[k] * xcoords[k + 1] ** 2 - a3 * xcoords[
                k] ** 2 * xcoords[k + 1]
            nu1 = a0 + a1 + a2 * (xcoords[k + 1] ** 2 + 2 * xcoords[k + 1] * xcoords[k]) + a3 * (
                    xcoords[k] ** 2 + 2 * xcoords[k] * xcoords[k + 1])
            nu2 = -a2 * (xcoords[k] + 2 * xcoords[k + 1]) - a3 * (xcoords[k + 1] + 2 * xcoords[k])
            nu3 = a2 + a3
            P = [nu0, nu1, nu2, nu3]
        end1 = time.time()
        time1 += [end1 - start1]
        # method 2
        start2 = time.time()
        xcoords = [a + k * (b - a) / n for k in range(n + 1)]
        ycoords = [f(k) for k in xcoords]
        dycoords = [misc.derivative(f, k) for k in xcoords]
        for k in range(n):
            M = np.array([(1, xcoords[k], xcoords[k] ** 2, xcoords[k] ** 3),
                          (1, xcoords[k + 1], xcoords[k + 1] ** 2, xcoords[k + 1] ** 3),
                          (0, 1, 2 * xcoords[k], 3 * xcoords[k] ** 2),
                          (0, 1, 2 * xcoords[k + 1], 3 * xcoords[k + 1] ** 2)])
            N = np.array([(ycoords[k]), (ycoords[k + 1]), (dycoords[k]), (dycoords[k + 1])])
            O = np.dot(np.linalg.inv(M), N)
            a0, a1, a2, a3 = O[0], O[1], O[2], O[3]
            P = [a0, a1, a2, a3]
        end2 = time.time()
        time2 += [end2 - start2]
        plt.clf()
    sigma1, sigma2 = 0, 0
    n = len(time1)
    print(time1, time2)
    for k in range(n):
        sigma1 += time1[k]
        sigma2 += time2[k]
    print(sigma1 / n, sigma2 / n)

# First method is a bit faster on reference functions
