from cubic_splines import interp_splines
from lagrange import tchebychev_integral
from polynomes import *

n = int(input("input n: "))
a = float(input("input a: "))
b = float(input("input b: "))

# Usual n-truncated power series
# exponential: the integral's value is exp(1)-1
exp = [1 / factorial(k) for k in range(n)]
I1 = polynomial_integral(exp, a, b)
I1_ref = np.exp(b) - np.exp(a)
# ch: the integral's value is sh(1)
ch = []
for k in range(n):
    if k % 2 == 0:
        ch.append(1 / factorial(k))
    else:
        ch.append(0)
I2 = polynomial_integral(ch, a, b)
I2_ref = np.sinh(b) - np.sinh(a)
# sh: the integral's value is ch(1)-1
sh = []
for k in range(n):
    if k % 2 == 1:
        sh.append(1 / factorial(k))
    else:
        sh.append(0)
I3 = polynomial_integral(sh, a, b)
I3_ref = np.cosh(b) - np.cosh(a)
# cos: the integral's value is sin(1)
cos = []
for k in range(n):
    if k % 2 == 0:
        cos.append((-1) ** (k / 2) / factorial(k))
    else:
        cos.append(0)
I4 = polynomial_integral(cos, a, b)
I4_ref = np.sin(b) - np.sin(a)
# sin: the integral's value is 1-cos(1)
sin = []
for k in range(n):
    if k % 2 == 1:
        sin.append((-1) ** ((k - 1) / 2) / factorial(k))
    else:
        sin.append(0)
I5 = polynomial_integral(sin, a, b)
I5_ref = np.cos(a) - np.cos(b)

# Tchebychev's n-degree method :
I1T = tchebychev_integral(np.exp, a, b, n)
I2T = tchebychev_integral(np.cosh, a, b, n)
I3T = tchebychev_integral(np.sinh, a, b, n)
I4T = tchebychev_integral(np.cos, a, b, n)
I5T = tchebychev_integral(np.sin, a, b, n)

# Cubic splines n-degree method:
I1S = interp_splines(np.exp, a, b, n)
I2S = interp_splines(np.cosh, a, b, n)
I3S = interp_splines(np.sinh, a, b, n)
I4S = interp_splines(np.cos, a, b, n)
I5S = interp_splines(np.sin, a, b, n)

print("Error SE :")
print(" -exp: ", abs(I1_ref - I1))
print(" -ch: ", abs(I2_ref - I2))
print(" -sh: ", abs(I3_ref - I3))
print(" -cos: ", abs(I4_ref - I4))
print(" -sin: ", abs(I5_ref - I5))
print("Error Tchebychev :")
print(" -exp: ", abs(I1_ref - I1T))
print(" -ch: ", abs(I2_ref - I2T))
print(" -sh: ", abs(I3_ref - I3T))
print(" -cos: ", abs(I4_ref - I4T))
print(" -sin: ", abs(I5_ref - I5T))
print("Error splines :")
print(" -exp: ", abs(I1_ref - I1S))
print(" -ch: ", abs(I2_ref - I2S))
print(" -sh: ", abs(I3_ref - I3S))
print(" -cos: ", abs(I4_ref - I4S))
print(" -sin: ", abs(I5_ref - I5S))
