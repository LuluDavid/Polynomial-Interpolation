# Polynomial_Interpolation

This git repository is a small implementation of polynomials and their use for interpolation.

* First of all, we define polynomials as lists of coefficients, and a few useful functions on them.
* Then, we define Lagrange interpolation for a given function.
* Then, we define Hermite interpolation for a given function.
* We also define Lagrange interpolation via Tchebychev points, in order to minimize the error between
  the interpolation Lagrange polynomial and the function it interpolates.
* Then, we introduce truncated power series and cubic splines to approximate functions (and integrals)
* Finally, we compare all those approaches with calculus on trivial functions (exp, sin, etc.)

This repo also contains a few attempts on those objects, such as 2D interpolation of colors on a picture.

The algorithms are quite optimized, even though the imports between different files are not resolved yet,
because I had no idea how to handle that back then.

My main usage for it was to successively compile every file in Pyzo, then call functions on different
functions.
