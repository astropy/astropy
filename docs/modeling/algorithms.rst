**********
Algorithms
**********

Univariate polynomial evaluation
================================

* The evaluation of 1-D polynomials uses Horner's algorithm.

* The evaluation of 1-D Chebyshev and Legendre polynomials uses Clenshaw's
  algorithm.


Multivariate polynomial evaluation
==================================

* Multivariate Polynomials are evaluated following the algorithm in [1]_ .  The
  algorithm uses the following notation:

  - **multiindex** is a tuple of non-negative integers for which the length is
    defined in the following way:

    .. math:: \alpha = (\alpha1, \alpha2, \alpha3),  |\alpha| = \alpha1+\alpha2+\alpha3


  - **inverse lexical order** is the ordering of monomials in such a way that
    :math:`{x^a < x^b}` if and only if there exists :math:`{1 \le i \le n}`
    such that :math:`{a_n = b_n, \dots, a_{i+1} = b_{i+1}, a_i < b_i}`.

    In this ordering :math:`y^2 > x^2*y` and :math:`x*y > y`

  - **Multivariate Horner scheme** uses d+1 variables :math:`r_0, ...,r_d` to
    store intermediate results, where *d* denotes the number of variables.

    Algorithm:

    1. Set *di* to the max number of variables (2 for a 2-D polynomials).

    2. Set :math:`r_0` to :math:`c_{\alpha(0)}`, where c is a list of
       coefficients for each multiindex in inverse lexical order.

    3. For each monomial, n, in the polynomial:

       - determine :math:`k = max \{1 \leq j \leq di: \alpha(n)_j \neq \alpha(n-1)_j\}`

       - Set :math:`r_k := l_k(x)* (r_0 + r_1 + \dots + r_k)`

       - Set :math:`r_0 = c_{\alpha(n)}, r_1 = \dots r_{k-1} = 0.`

    4. return :math:`r_0 + \dots + r_{di}`

* The evaluation of multivariate Chebyshev and Legendre polynomials uses a
  variation of the above Horner's scheme, in which every Legendre or Chebyshev
  function is considered a separate variable.  In this case the length of the
  :math:`\alpha` indices tuple is equal to the number of functions in x plus
  the number of functions in y.  In addition the Chebyshev and Legendre
  functions are cached for efficiency.



.. [1] J. M. Pena, Thomas Sauer, "On the Multivariate Horner Scheme", SIAM Journal on Numerical Analysis, Vol 37, No. 4
