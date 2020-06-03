# Lagrange interpolation, using Tchebychev abscisses, Newton-Côtes formulas

from polynomes import *


# Compute Li
# n**3 operations
def lagrange_polynomial_ith(xcoords, i):
    li = [1]
    n = len(xcoords)
    for k in range(0, n):
        if k != i:
            Pj = [-xcoords[k] / (xcoords[i] - xcoords[k]), 1 / (xcoords[i] - xcoords[k])]
            li = multiply_polynomials(li, Pj)
    return li


# Compute the complete lagrange polynomial
# n**4 operations
def lagrange_polynomial(xcoords, ycoords):
    L = [0]
    n = len(xcoords)
    for i in range(0, n):
        Li = lagrange_polynomial_ith(xcoords, i)
        Li = multiply_polynomial(Li, ycoords[i])
        L = add_polynomials(L, Li)
    return L


# The polynomial's graph
def lagrange_graph(xcoords, ycoords):
    L = lagrange_polynomial(xcoords, ycoords)
    a, b = min(xcoords), max(xcoords)
    coordsx = np.linspace(a, b, 1000)
    coordsy = [polynomial_image(L, x) for x in coordsx]
    plt.plot(coordsx, coordsy, "r-")
    plt.show()


# The polynomial's graph for ycoords = map f xcoords
def lagrange_graph_function(xcoords, f):
    ycoords = [f(x) for x in xcoords]
    lagrange_graph(xcoords, ycoords)


# Compute Lagrange's n degree polynomial with Tchebychev's abscisses
# n**4 operations
def lagrange_tchebychev(f, a, b, n):
    n += 1
    racines = [(a + b) / 2 + (a - b) / 2 * np.cos((2 * k + 1) * np.pi / (2 * n)) for k in range(0, n)]
    ycoords = [f(x) for x in racines]
    L = lagrange_polynomial(racines, ycoords)
    return L


def lagrange_tchebychev_graph(f, a, b, n):
    L = lagrange_tchebychev(f, a, b, n)
    coordsx = np.linspace(a, b, 1000)
    coords = [f(x) for x in coordsx]
    coordsy = [polynomial_image(L, x) for x in coordsx]
    plt.clf()
    plt.plot(coordsx, coords, "b-", linewidth=2)
    plt.plot(coordsx, coordsy, "r-")
    plt.legend(loc="best")
    plt.show()


# Comparing a random n degree Lagrange polynomial with Tchebychev
def random_vs_tchebychev(f, n, a, b):
    res = []
    for k in range(0, n + 1):
        res += [np.random.uniform(a, b)]
    ycoords = [f(x) for x in res]
    L = lagrange_polynomial(res, ycoords)
    coordsx = np.linspace(a, b, 1000)
    coords = [f(x) for x in coordsx]
    coordsy = [polynomial_image(L, x) for x in coordsx]
    plt.clf()
    plt.plot(coordsx, coords, "b-", linewidth=2)
    plt.plot(coordsx, coordsy, "r--", label="lagrange's interpolation polynomial")
    lagrange_tchebychev_graph(f, a, b, n)


# Newton-Côtes formula using Lagrange's polynomials for integral interpolation
# n**4 operations
def newton_cotes(f, a, b, n):
    xcoords = [a + k * (b - a) / n for k in range(0, n)]
    sigma = 0
    for i in range(0, n):
        li = lagrange_polynomial_ith(xcoords, i)
        sigma += f(xcoords[i]) * polynomial_integral(li, a, b)
    return sigma


# Approximation of the function's integral with Tchebychev's n degree polynomial
# n**4 operations
def tchebychev_integral(f, a, b, n):
    L = lagrange_tchebychev(f, a, b, n)
    return polynomial_integral(L, a, b)
