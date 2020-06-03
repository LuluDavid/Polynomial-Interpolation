# Polynomial_Interpolation

This git repository is a small implementation of polynomials and their use for interpolation, in the context
of a work for undergraduate courses in France (TIPE in Prepa).

The final report (in french) is accessible [here](https://drive.google.com/file/d/1UaUzbZDAs9JtF5YragzRGSl5GUYIQsT2/view?usp=sharing), I might do a latex version on day.

## Polynomials

I introduced polynomials as lists of coefficients, and defined on them in **polynomials.py**. A few utils are also introduced there : polynomial sum, polynomial constant, polynomial product, polynomial image, polynomial integral and derivate function, polynomial display with *matplotlib.pyplot*.

## Lagrange

In the file **Lagrange.py**, we introduce the Lagrange polynomials, and their representation. Here are the function's specs :
  * *lagrange_polynomial_i(xcoords, i)* : calculate Lagrange's ith polynomial in the Lagrange's formula for xcoords
  * *lagrange_polynomial(xcoords, ycoords)* : calculate Lagrange's polynomial for xcoords, ycoords
  * *lagrange_graph(xcoords, ycoords)* : calculate and then display Lagrange's polynomial on xcoords
  * *lagrange_graph_function(xcoords, f)* : do the preveious operation with ```ycoords = map f xcoords```
  * *lagrange_tchebychev(f, a, b, n)* : calculate Tchebychev's coordinates (for degree n), and calculate Lagrange's polynomial for these xcoords (they should limit the polynomial's sup between a and b for all n degree polynomials)
  * *lagrange_tchebychev_graph(f, a, b, n)* : display the previous polynomial on generated xcoords
  * *newton_cotes(f, a, b, n)* : calculate the n-th Newton-Côtes formula thanks to Lagrange's polynomials (integral approximation)
  * *tchebychev_integral(f, a, b, n)* : calculate the function's integral approximation thanks to its n-th degree tchebychev polynomial calculated between a and b

## L'Hermite

In the file **hermite.py**, it is quite similar to Lagrange's functions, except this time it is done for Hermite's interpolation

## Points interpolation

In the file **scatter_interpolation.py**, I defined a simple function *lagrange_tchebychev_pts(xcoords, ycoords, n)* to calculate the nearest points to n degree Tchebychev points between min and max of xcoords, and interpolate those points with a Lagrange polynomial, to then display it. 

The goal is to define a good polynomial interpolation for a set of points.

## 2D interpolation

In the file **lagrange2D.py**, I tried to simply interpolation according to two directions, x and y, and then applied it to a picture where I set some colors at the extremal points of the picture, and then interpolated the colors between them, which gave the following result.

I could try it with l'Hermite's interpolation one day.

## Cubic Splines

In the file **cubic_splines.py**, I defined *cubic_splines(f, a, b, n)*, which displays the n subdivision cubic splines interpolation on function f between a and b (and an alternative definition *cubic_splines2*, as well as in another file, but it is not that important).

I also defined a function *interp_splines(f, a, b, n)* to approximate the integral of f between a and b with the sum of the splines integrals.

## Comparison

I compared all those methods, plus *truncated power series*, in the file **integral_comparison.py**, on well known functions, such as the exponential, cosh, sinh, cos and sin.

You just need to run the file (after having run the other files to access their definitions), and then choose the degree n, and the interval [a , b].

## TODO

As you can see, this project is quite old, and there are no imports, because I had a minimal knowledge of python.
Moreover, it's written entirely in french, which is not a good thing.
Finally, the associated report is also in french, which does not help to understand.
I might work at it when I have time, but for the moment it will stay like that.

If you have any recommandation, do not hesitate to message me or to submit a pull request.
