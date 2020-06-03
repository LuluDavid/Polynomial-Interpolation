import numpy as np
import matplotlib.pyplot as plt


# Interpolating missing pixels from the 4 extremal colors
def bilinear_interpolation(mat, m):
    l, p = np.shape(mat)
    result = np.zeros((l * (m + 1) - m, p * (m + 1) - m))
    for u in range(l):
        for v in range(p):
            result[(m + 1) * u][(m + 1) * v] = mat[u][v]
    k0, l0 = 0, 0
    a, b, c, d = 0, 0, 0, 0
    while k0 <= l - 2:
        l0 = 0
        while l0 <= p - 2:
            a, b, c, d = mat[k0][l0], mat[k0][l0 + 1], mat[k0 + 1][l0], mat[k0 + 1][l0 + 1]
            for u in range(k0, k0 + m + 1):
                for v in range(l0, l0 + m + 1):
                    if (u, v) != (k0, l0) and (u, v) != (k0, l0 + m) and (u, v) != (k0 + m, l0) and (u, v) != (
                            k0 + m, l0 + m):
                        result[u][v] = f(u, v, a, b, c, d)
            l0 += 1
        for u in range(k0, k0 + m):
            if u != k0 and u != k0 + m:
                result[u][p * (m + 1) - m - 1] = f(u, p * (m + 1) - m - 1, a, b, c, d)
        k0 += 1
    for v in range(l0, l0 + m):
        if v != l0 and v != l0 + m:
            result[l * (m + 1) - m - 1][v] = f(l * (m + 1) - m - 1, v, a, b, c, d)
    return result


n = int(input("Which with for the square ?"))
res = [[[0, 0, 0] for k in range(n)] for l in range(n)]
for k in range(n):
    res[0][k], res[n - 1][k] = [100, 0, 0], [0, 0, 100]

f0, f1, f2, f3 = res[0][0], res[0][n - 1], res[n - 1][0], res[n - 1][n - 1]


def f(i, j, a, b, c, d):
    result = ((a * (i - (n - 1)) - c * i) * (j - (n - 1)) + (b * (i - (n - 1)) - d * i) * j) / ((n - 1) ** 2)
    return result


w0, x0, y0, z0 = f0[0], f1[0], f2[0], f3[0]
w1, x1, y1, z1 = f0[1], f1[1], f2[1], f3[1]
w2, x2, y2, z2 = f0[2], f1[2], f2[2], f3[2]
for i in range(1, n - 1):
    for j in range(n):
        if (i, j) != (0, 0) and (i, j) != (0, n - 1) and (i, j) != (n - 1, 0) and (i, j) != (n - 1, n - 1):
            res[i][j] = [f(i, j, w0, x0, y0, z0), f(i, j, w1, x1, y1, z1), f(i, j, w2, x2, y2, z2)]
plt.imshow(res)
plt.show()
