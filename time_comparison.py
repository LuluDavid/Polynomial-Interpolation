import time

import numpy as np
from scipy import misc

from hermite import hermite_polynomial
from lagrange import lagrange_polynomial
from power_series import power_serie

f = np.sin
a = -5
b = 5
n_ = int(input("input n ? (polynomial's degree will be 2n+1)"))
startLT = time.time()
n = 2 * n_ + 2
racines = [(a + b) / 2 + (a - b) / 2 * np.cos((2 * k + 1) * (np.pi) / (2 * n)) for k in range(0, n)]
ycoords = [f(x) for x in racines]
L = lagrange_polynomial(racines, ycoords)
endLT = time.time()
lengthLT = endLT - startLT
startLH = time.time()
n_ = n_ + 1
ycoords = [f(x) for x in racines]
dycoords = [misc.derivative(f, x) for x in racines]
L = hermite_polynomial(racines, ycoords, dycoords)
endLH = time.time()
lengthLH = endLH - startLH
startSP = time.time()
res = []
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
    res.append(P)
endSP = time.time()
lengthSP = endSP - startSP
startSE = time.time()
L = power_serie(f, n)
endSE = time.time()
lengthSE = endSE - startSE
print("-> Lagrange's length:", lengthLH)
print("-> L'Hermite's length:", lengthLT)
print("-> Cubic splines's length:", lengthSP)
print("-> Power series's length", lengthSE)
