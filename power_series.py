from decimal import *
import numpy as np

from polynomes import factorial


def derivativen(f, x, n, epsilon=0.001):
    M = [[k ** j for k in range(n + 1)] for j in range(n + 1)]
    B = [[0] for k in range(n)] + [[1]]
    res = system_solving(M, B)
    ftot = Decimal(0)
    for k in range(n + 1):
        fk = Decimal(f(x + k * epsilon))
        ftot += Decimal(res[k][0]) * fk
    ftot = ftot * Decimal(factorial(n)) / (Decimal(epsilon) ** n)
    return ftot


def power_serie(f, n):
    P = [derivativen(f, 0, k) for k in range(n)]
    return P


# Matrix inversion (for Vandermonde)
def lines_permutation(A, i, j):
    p = np.shape(A)[1]
    for k in range(p):
        A[i][k], A[j][k] = A[j][k], A[i][k]


def system_solving(A0, B0):
    A = np.copy(A0)
    Y = np.copy(B0)
    n, p = np.shape(A)
    nb_1, nb_c = np.shape(Y)
    for j in range(n - 1):
        m = abs(A[j, j])
        i_pivot = j
        for i in range(j + 1, n):
            if abs(A[i, j]) > m:
                m = abs(A[i, j])
                i_pivot = i
        lines_permutation(A, i_pivot, j)
        lines_permutation(Y, i_pivot, j)
        if A[j, j] != 0:
            for k in range(j + 1, n):
                Y[k, :] -= A[k, j] / A[j, j] * Y[j, :]
                A[k, :] -= A[k, j] / A[j, j] * A[j, :]
    X = np.zeros((n, nb_c))
    for i in range(n - 1, -1, -1):
        for k in range(i + 1, n):
            Y[i, :] -= A[i, k] * X[k, :]
        X[i, :] = Y[i, :] / A[i, i]
    return X
