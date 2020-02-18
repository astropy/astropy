.. include:: links.inc
	     
.. _polynomial_models:

*****************
Polynomial Models
*****************

.. _domain-window-note:

Notes regarding usage of domain and window
------------------------------------------

Most of the polynomial models have optional domain and window attributes.
It is important to understand how they currently are interpreted, which
can be confusing since the terminology often implies something different.

Both the domain and window attributes for a polynomial consist of a two
element list (this will change to tuples in a future release) that
indicate a range of values for input values. For 2-Dimensional polynomials
the attributes become x_domain, y_domain, x_window, and y_window.
Generally speaking, the main purpose of these attributes is to define
a linear transform between the supplied input variable and the resultant
input variable that is supplied to the polynomial. For example, if
domain = [-2, 2] and window = [-1, 1], input values will be divided by
two so that the domain maps to the window. Correspondingly the pair
domain = [0, 2], window = [-1, 1] implies that 1 will be subtracted from
the input variable before using it in the polynomial.

Neither domain or window are meant to imply that values that fall outside
of their corresponding ranges will result in an exception, or that
such values are necessarily invalid (the latter depends on the context
of how the polynomial is being used).

It is the case that the orthogonal polynomials are defined on a range of
[-1, 1], but nothing in the current machinery prevents them from being
evaluated outside that range.

Domain is used in fitting polynomials to bound the input variable to map
to the defined window so that they fall within the expected [-1, 1] range
for such polynomials. That is, the fitting routine will set the domain to
map to the window range for the range of input x values supplied (so that
domain may change if the minimum and maximum x values being fit change).

The meaning of these terms may conflict with expectations (e.g., domain
is often meant to mean the range of input values the funciton is valid
for). For fit results that is somewhat true, but otherwise, it is not.
The default values for ordinary polynomials is [-1, 1] for both domain
and window, which effectively signals no transformation of the input
variable.

The terminology was adopted from numpy polynomials, which have the same
confusion in meaning.


**************
1D Polynomials
**************

- :class:`~astropy.modeling.polynomial.Polynomial1D`

- :class:`~astropy.modeling.polynomial.Chebyshev1D`

- :class:`~astropy.modeling.polynomial.Legendre1D`

- :class:`~astropy.modeling.polynomial.Hermite1D`

**************
2D Polynomials
**************

- :class:`~astropy.modeling.polynomial.Polynomial2D`

- :class:`~astropy.modeling.polynomial.Chebyshev2D`

- :class:`~astropy.modeling.polynomial.Legendre2D`

- :class:`~astropy.modeling.polynomial.Hermite2D`

- :class:`~astropy.modeling.polynomial.SIP` model implements the
  Simple Imaging Polynomial (`SIP`_) convention 

  
