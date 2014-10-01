****************************
Image stretching and scaling
****************************

The `imageutils.scaling` module provides a framework for dealing with image
stretching and scaling, as well as determining intervals based on various
criteria (such as percentile values).

Intervals
=========

Several classes are provided for determining intervals and for scaling values
in this interval to the [0:1] range. One of the simplest examples is the
:class:`~imageutils.scaling.MinMaxInterval` which determines the limits of the
values based on the minimum and maximum values in the array. The class is
instantiated with no arguments::


    >>> from imageutils.scaling import MinMaxInterval
    >>> interval = MinMaxInterval()

and the limits can be determined by calling the
:meth:`~imageutils.scaling.MinMaxInterval.get_limits` method, which takes the
array of values::

    >>> vmin, vmax = interval.get_limits([1, 3, 4, 5, 6])
    >>> vmin
    1
    >>> vmax
    6

The ``interval`` instance can also be called like a function to actually
normalize values to the range::

    >>> interval([1, 3, 4, 5, 6])
    array([ 0. ,  0.4,  0.6,  0.8,  1. ])

Other interval classes include :class:`~imageutils.scaling.ManualInterval`,
:class:`~imageutils.scaling.PercentileInterval`, and
:class:`~imageutils.scaling.AsymmetricPercentileInterval`. For these three,
values in the array can fall outside of the limits given by the interval. A
``clip`` argument is provided to control the behavior of the scaling when
values fall outside the limits::


    >>> from imageutils.scaling import PercentileInterval
    >>> interval = PercentileInterval(50.)
    >>> vmin, vmax = interval.get_limits([1, 3, 4, 5, 6])
    >>> vmin
    3.0
    >>> vmax
    5.0
    >>> interval([1, 3, 4, 5, 6])  # default is clip=True
    array([ 0. ,  0. ,  0.5,  1. ,  1. ])
    >>> interval([1, 3, 4, 5, 6], clip=False)
    array([-1. ,  0. ,  0.5,  1. ,  1.5])

Stretching
==========

In addition to classes that can scale values to the [0:1] range, a number of
classes are provide to 'stretch' the values using different functions. These
map a [0:1] range onto a transformed [0:1] range. A simple example is the
:class:`~imageutils.scaling.SqrtStretch` class::

    >>> from imageutils.scaling import SqrtStretch
    >>> stretch = SqrtStretch()
    >>> stretch([0., 0.25, 0.5, 0.75, 1.])
    array([ 0.        ,  0.5       ,  0.70710678,  0.8660254 ,  1.        ])

As for the intervals, values outside the [0:1] range can be treated differently
depending on the ``clip`` argument. By default, output values are clipped to
the [0:1] range::


    >>> stretch([-1., 0., 0.5, 1., 1.5])
    array([ 0.       ,  0.        ,  0.70710678,  1.        ,  1.        ])

but this can be disabled::

    >>> stretch([-1., 0., 0.5, 1., 1.5], clip=False)
    array([        nan,  0.        ,  0.70710678,  1.        ,  1.22474487])

.. note:: The stretch functions are similar but not always strictly identical
          to those used in e.g. ds9 (although they should have the same
          behavior). The equations for the ds9 stretches can be found `here
          <http://ds9.si.edu/ref/how.html>`_ and can be compared to the
          equations for our stretches provided in the `imageutils.scaling` API
          section. The main difference between our stretches and ds9 is that we
          have adjusted them so that the [0:1] range always maps exactly to the
          [0:1] range.

Combining transformations
=========================

Any stretches and intervals can be chained by using the ``+`` operator, which
returns a new transformation. For example, to apply scaling based on a
percentile value, followed by a square root stretch, you can do::

    >>> transform = SqrtStretch() + PercentileInterval(90.)
    >>> transform([1, 3, 4, 5, 6])
    array([ 0.        ,  0.60302269,  0.76870611,  0.90453403,  1.        ])
    
As before, the combined transformation can also accept a ``clip`` argument
(which is `True` by default).

.. automodapi:: imageutils.scaling
