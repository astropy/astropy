.. _quantity:

Plotting Astropy objects in Matplotlib
**************************************

.. |quantity| replace:: :class:`~astropy.units.Quantity`
.. |time| replace:: :class:`~astropy.time.Time`

.. _plotting-quantities:

Plotting quantities
===================

|quantity| objects can be conveniently plotted using matplotlib.  This
feature needs to be explicitly turned on:

.. doctest-requires:: matplotlib

    >>> from astropy.visualization import quantity_support
    >>> quantity_support()  # doctest: +IGNORE_OUTPUT
    <astropy.visualization.units.MplQuantityConverter ...>

Then |quantity| objects can be passed to matplotlib plotting
functions.  The axis labels are automatically labeled with the unit of
the quantity:

.. doctest-requires:: matplotlib

    >>> from matplotlib import pyplot as plt
    >>> plt.figure(figsize=(5,3))
    <...>
    >>> plt.plot([1, 2, 3] * u.m)
    [...]

.. plot::

    from astropy import units as u
    from astropy.visualization import quantity_support
    quantity_support()
    from matplotlib import pyplot as plt
    plt.figure(figsize=(5,3))
    plt.plot([1, 2, 3] * u.m)

Quantities are automatically converted to the first unit set on a
particular axis, so in the following, the y-axis remains in ``m`` even
though the second line is given in ``cm``:

.. doctest-requires:: matplotlib

    >>> plt.plot([1, 2, 3] * u.cm)
    [...]

.. plot::

    from astropy import units as u
    from astropy.visualization import quantity_support
    quantity_support()
    from matplotlib import pyplot as plt
    plt.figure(figsize=(5,3))
    plt.plot([1, 2, 3] * u.m)
    plt.plot([1, 2, 3] * u.cm)

Plotting a quantity with an incompatible unit will raise an exception.
For example, calling ``plt.plot([1, 2, 3] * u.kg)`` (mass unit) to overplot
on the plot above that is displaying length units.

To make sure unit support is turned off afterward, you can use
`~astropy.visualization.quantity_support` with a ``with`` statement:

.. doctest-requires:: matplotlib

    >>> from astropy.visualization import quantity_support
    >>> from matplotlib import pyplot as plt
    >>> with quantity_support():
    ...     plt.figure(figsize=(5,3))
    ...     plt.plot([1, 2, 3] * u.m)
    <...>
    [...]

.. plot::

    from astropy import units as u
    from astropy.visualization import quantity_support
    from matplotlib import pyplot as plt
    with quantity_support():
        plt.figure(figsize=(5,3))
        plt.plot([1, 2, 3] * u.m)

.. _plotting-times:

Plotting times
==============

Similarly to |quantity|, |time| objects can also be plotted using matplotlib
in a way that the scale and format used for the axes can be controlled. This
feature needs to be explicitly turned on:

.. doctest-requires:: matplotlib

    >>> from astropy.visualization import time_support
    >>> time_support()  # doctest: +IGNORE_OUTPUT
    <astropy.visualization.units.MplTimeConverter ...>

Once this is enabled, |time| objects can be passed to matplotlib plotting
functions. The axis labels are then automatically labeled with times formatted
using the |time| class:

.. doctest-requires:: matplotlib

    >>> from matplotlib import pyplot as plt
    >>> from astropy.time import Time
    >>> plt.figure(figsize=(5,3))
    <...>
    >>> plt.plot(Time(...))
    [...]

.. plot::

  :include-source:
  :context: reset

    from matplotlib import pyplot as plt
    from astropy.time import Time
    from astropy.visualization import time_support

    time_support()

    plt.figure(figsize=(5,3))
    plt.plot(Time([58000, 59000, 62000], format='mjd'))

By default, times are shown in UTC and in the ISO format, but this can be
controlled by passing arguments to ``time_support``::

  .. plot::

    :include-source:
    :context:

      time_support(format='mjd', scale='tai')
      plt.figure(figsize=(5,3))
      plt.plot(Time([58000, 59000, 62000], format='mjd'))

To make sure support for plotting times is turned off afterward, you can use
`~astropy.visualization.time_support` as a context manager::

    with time_support(format='mjd', scale='tai'):
        plt.figure(figsize=(5,3))
        plt.plot(Time([58000, 59000, 62000], format='mjd'))
