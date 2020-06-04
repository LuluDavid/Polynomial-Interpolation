# Basic functions on polynomials

import numpy as np
import matplotlib.pyplot as plt

# Polynomials as lists of coefficients

# P+Q: max(deg(P),deg(Q)) operations


def add_polynomials(P, Q):
    somme = []
    n = len(P)
    m = len(Q)
    if n < m:
        P += (m - n) * [0]
    elif n < m:
        Q += (n - m) * [0]
    deg_max = max(n, m)
    for k in range(0, deg_max):
        somme += [P[k] + Q[k]]
    return somme


# aP: deg(P) operations
def multiply_polynomial(P, a):
    n = len(P)
    for k in range(0, n):
        P[k] *= a
    return P


# PQ: (deg(P)+deg(Q))**2 operations
def multiply_polynomials(P, Q):
    n = len(P)
    m = len(Q)
    P += m * [0]
    Q += n * [0]
    p = n + m
    # thus deg(PQ)=n+m
    res = []
    for k in range(0, p - 1):
        a_k = 0
        for i in range(0, p - 1):
            a_k += P[i] * Q[k - i]
        res += [a_k]
    return res


# P(x): deg(P) operations
def polynomial_image(P, x):
    n = len(P)
    sigma = 0
    for k in range(0, n):
        sigma += P[k] * x ** k
    return sigma


# int(P): 3deg(P) operations
def polynomial_integral(P, a, b):
    n = len(P)
    Q = []
    for k in range(0, n):
        Q += [P[k] / (k + 1)]
    Q = [0] + Q
    I = polynomial_image(Q, b) - polynomial_image(Q, a)
    return I


# P: deg(P) operations
def polynomial_derivative(P):
    n = len(P)
    dP = []
    for k in range(1, n):
        dP += [P[k] * (k)]
    return dP


def polynomial_graph(P, a, b):
    coordsx = np.linspace(a, b, 1000)
    coordsy = [polynomial_image(P, x) for x in coordsx]
    plt.clf()
    plt.plot(coordsx, coordsy, "r-", linewidth=2)


# Two useful utils to compute power series : factorial and binomial coeffs


# log(n) operations
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n - 1)


# 3*log(n) operations
def binomial_coefficients(k, n):
    u = factorial(n) / (factorial(k) * factorial(n - k))
    return u
