.. _structured_units:

Structured Units
****************

Numpy arrays can be :doc:`structured arrays <numpy:user/basics.rec>`, where
each element consists of multiple fields. These can be used with |Quantity|
using a |StructuredUnit|, which provides a |Unit| for each field. For example,
this allows constructing a single |Quantity| object with position and velocity
fields that have different units, but are contained within the same object
(as is needed to support units in the |PyERFA| wrappers around the |ERFA|
routines that use position-velocity arrays).

Creating Structured Quantities
==============================

You can create structured quantities either directly or by multiplication with
a |StructuredUnit|, with the latter in turn either created directly, or
through `~astropy.units.Unit`.

Example
-------

.. EXAMPLE START: Creating Structured Quantities

To create a structured quantity containing a position and velocity::

  >>> import astropy.units as u, numpy as np
  >>> pv_values = np.array([([1., 0., 0.], [0., 0.125, 0.]),
  ...                       ([0., 1., 0.], [-0.125, 0., 0.])],
  ...                      dtype=[('p', '(3,)f8'), ('v', '(3,)f8')])
  >>> pv = u.Quantity(pv_values, u.StructuredUnit((u.km, u.km/u.s)))
  >>> pv
  <Quantity [([1., 0., 0.], [ 0.   ,  0.125,  0.   ]),
             ([0., 1., 0.], [-0.125,  0.   ,  0.   ])] (km, km / s)>
  >>> pv_values * u.Unit('AU, AU/day')
  <Quantity [([1., 0., 0.], [ 0.   ,  0.125,  0.   ]),
             ([0., 1., 0.], [-0.125,  0.   ,  0.   ])] (AU, AU / d)>

As for normal |Quantity| objects, you can access the value and the unit with the
`~astropy.units.Quantity.value` and `~astropy.units.Quantity.unit` attribute,
respectively. In addition, you can index any given field using its name::

  >>> pv = pv_values * u.Unit('km, km/s')
  >>> pv.value
  array([([1., 0., 0.], [ 0.   ,  0.125,  0.   ]),
         ([0., 1., 0.], [-0.125,  0.   ,  0.   ])],
        dtype=[('p', '<f8', (3,)), ('v', '<f8', (3,))])
  >>> pv.unit
  Unit("(km, km / s)")
  >>> pv['v']
  <Quantity [[ 0.   ,  0.125,  0.   ],
             [-0.125,  0.   ,  0.   ]] km / s>

Structures can be nested, as in this example taken from an |PyERFA| test case
for :func:`erfa.ldn`::

  >>> ldbody = [
  ...     (0.00028574, 3e-10, ([-7.81014427, -5.60956681, -1.98079819],
  ...                          [0.0030723249, -0.00406995477, -0.00181335842])),
  ...     (0.00095435, 3e-9, ([0.738098796, 4.63658692, 1.9693136],
  ...                         [-0.00755816922, 0.00126913722, 0.000727999001])),
  ...     (1.0, 6e-6, ([-0.000712174377, -0.00230478303, -0.00105865966],
  ...                  [6.29235213e-6, -3.30888387e-7, -2.96486623e-7]))
  ...     ] * u.Unit('Msun,radian,(AU,AU/day)')
  >>> ldbody  # doctest: +FLOAT_CMP
  <Quantity [(2.8574e-04, 3.e-10, ([-7.81014427e+00, -5.60956681e+00, -1.98079819e+00], [ 3.07232490e-03, -4.06995477e-03, -1.81335842e-03])),
             (9.5435e-04, 3.e-09, ([ 7.38098796e-01,  4.63658692e+00,  1.96931360e+00], [-7.55816922e-03,  1.26913722e-03,  7.27999001e-04])),
             (1.0000e+00, 6.e-06, ([-7.12174377e-04, -2.30478303e-03, -1.05865966e-03], [ 6.29235213e-06, -3.30888387e-07, -2.96486623e-07]))] (solMass, rad, (AU, AU / d))>

.. EXAMPLE END

Converting to Different Units
=============================

Like regular |Quantity| objects, structured quantities can be converted to
different units, as long as they have the same structure and each unit is
equivalent.

Example
-------

.. EXAMPLE START: Converting Structured Quantities to Different Units

To convert a structured quantity to a different unit::

  >>> pv.to((u.m, u.m / u.s))  # doctest: +FLOAT_CMP
  <Quantity [([1000.,    0.,    0.], [   0.,  125.,    0.]),
             ([   0., 1000.,    0.], [-125.,    0.,    0.])] (m, m / s)>
  >>> pv.cgs
  <Quantity [([100000.,      0.,      0.], [     0.,  12500.,      0.]),
             ([     0., 100000.,      0.], [-12500.,      0.,      0.])] (cm, cm / s)>

.. EXAMPLE END

Use with ERFA
=============

The |ERFA| C routines make use of structured types, and these are exposed in
the |PyERFA| interface.

.. warning:: Not all |PyERFA| routines are wrapped yet. Help with adding
             wrappers will be appreciated.

Example
-------

.. EXAMPLE START: Using Structured Quantities with ERFA

To use a position-velocity structured array with |PyERFA|::

  >>> import erfa
  >>> pv_values = np.array([([1., 0., 0.], [0., 0.125, 0.]),
  ...                       ([0., 1., 0.], [-0.125, 0., 0.])],
  ...                      dtype=erfa.dt_pv)
  >>> pv = pv_values << u.Unit('AU,AU/day')
  >>> erfa.pvu(86400*u.s, pv)
  <Quantity [([ 1.   ,  0.125,  0.   ], [ 0.   ,  0.125,  0.   ]),
             ([-0.125,  1.   ,  0.   ], [-0.125,  0.   ,  0.   ])] (AU, AU / d)>
  >>> erfa.pv2s(pv)  # doctest: +FLOAT_CMP
  (<Quantity [0.        , 1.57079633] rad>,
   <Quantity [0., 0.] rad>,
   <Quantity [1., 1.] AU>,
   <Quantity [0.125, 0.125] rad / d>,
   <Quantity [0., 0.] rad / d>,
   <Quantity [0., 0.] AU / d>)
  >>> z_axis = np.array(([0, 0, 1], [0, 0, 0]), erfa.dt_pv) * u.Unit('1,1/s')
  >>> erfa.pvxpv(pv, z_axis)
  <Quantity [([ 0., -1.,  0.], [0.125, 0.   , 0.   ]),
             ([ 1.,  0.,  0.], [0.   , 0.125, 0.   ])] (AU, AU / d)>

.. EXAMPLE END
