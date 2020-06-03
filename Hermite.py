# Hermite interpolation


# Calculate hermite's interpolation polynomial from x, f(x) and f'(x) lists
# n**5 operations

from lagrange import *
from scipy import misc


def hermite_polynomial(xcoords, ycoords, dycoords):
    P = []
    n = len(xcoords)
    for i in range(n):
        Li = lagrange_polynomial_ith(xcoords, i)
        dLi = polynomial_derivative(Li)
        dLi_xi = polynomial_image(dLi, xcoords[i])
        hi = multiply_polynomials(add_polynomials([1], multiply_polynomial([2 * xcoords[i], -2], dLi_xi)),
                                  multiply_polynomials(Li, Li))
        ki = multiply_polynomials([-xcoords[i], 1], multiply_polynomials(Li, Li))
        P = add_polynomials(P,
                            add_polynomials(multiply_polynomial(hi, ycoords[i]), multiply_polynomial(ki, dycoords[i])))
    return P


def hermite_graph(xcoords, ycoords, dycoords):
    P = hermite_polynomial(xcoords, ycoords, dycoords)
    a, b = min(xcoords), max(xcoords)
    coordsx = np.linspace(a, b, 1000)
    coordsy = [polynomial_image(P, x) for x in coordsx]
    plt.plot(coordsx, coordsy, "r-")


def hermite_graph_function(xcoords, f):
    ycoords = [f(x) for x in xcoords]
    dycoords = [misc.derivative(f, x) for x in xcoords]
    hermite_graph(xcoords, ycoords, dycoords)


def hermite_tchebychev(f, a, b, n):
    n += 1
    racines = [(a + b) / 2 + (a - b) / 2 * np.cos((2 * k + 1) * (np.pi) / (2 * n)) for k in range(0, n)]
    ycoords = [f(x) for x in racines]
    dycoords = [misc.derivative(f, x) for x in racines]
    L = hermite_polynomial(racines, ycoords, dycoords)
    return L


def hermite_graph_tchebychev(f, a, b, n):
    L = hermite_tchebychev(f, a, b, n)
    coordsx = np.linspace(a, b, 1000)
    coords = [f(x) for x in coordsx]
    coordsy = [polynomial_image(L, x) for x in coordsx]
    plt.clf()
    plt.plot(coordsx, coords, "r-", linewidth=2, label="function to interpolate")
    plt.plot(coordsx, coordsy, "b-", label="tchebychev's n degree interpolating polynomial")
    plt.legend(loc="best")
    plt.show()


def random_vs_tchebychev_h(f, a, b, n):
    res = []
    LT = hermite_tchebychev(f, a, b, n)
    for k in range(0, n + 1):
        res += [np.random.uniform(a, b)]
    ycoords = [f(x) for x in res]
    dycoords = [misc.derivative(f, x) for x in res]
    L = hermite_polynomial(res, ycoords, dycoords)
    coordsx = np.linspace(a, b, 1000)
    coords = [f(x) for x in coordsx]
    coordsy = [polynomial_image(L, x) for x in coordsx]
    coordsy2 = [polynomial_image(LT, x) for x in coordsx]
    plt.clf()
    plt.plot(coordsx, coords, "b-", label="function to interpolate")
    plt.plot(coordsx, coordsy, "r--", label="hermite's interpolating polynomial")
    plt.plot(coordsx, coordsy2, "g--", label="hermite's tchebychev interpolating polynomial")
    plt.legend(loc="best")
    plt.show()


def tchebychev_integral_h(f, a, b, n):
    L = hermite_tchebychev(f, a, b, n)
    return polynomial_integral(L, a, b)
