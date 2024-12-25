from astropy.visualization.time import time_support
from astropy.visualization.units import quantity_support
from contextlib import contextmanager, ExitStack


@contextmanager
def astropy_support(*, quantity_support_kwargs=None, time_support_kwargs=None):
    """
    Enable support for plotting `astropy.units.Quantity` and `astropy.time.Time` instances in
    matplotlib.

    It can be used as a decorator or with a ``with`` statement.

    Examples
    ----------

    >>> import matplotlib.pyplot as plt
    >>> import astropy.units as u
    >>> from astropy.time import Time
    >>> from astropy.visualization.astropy_support import astropy_support

    >>> @astropy_support()
    ... def plot_example():
    ...     plt.figure()
    ...     plt.plot([1, 2, 3] * u.m)
    ...     plt.plot(Time(['2000-01-01', '2000-01-02', '2000-01-03']).plot_date)
    ...     plt.draw()
    ...     plt.show()

    >>> with astropy_support():  # doctest: +IGNORE_OUTPUT
    ...     plt.figure()
    ...     plt.plot([1, 2, 3] * u.m)
    ...     plt.plot(Time(['2000-01-01', '2000-01-02', '2000-01-03']).plot_date)
    ...     plt.draw()

    """
    with ExitStack() as stack:
        stack.enter_context(quantity_support(**(quantity_support_kwargs or {})))
        stack.enter_context(time_support(**(time_support_kwargs or {})))
        yield
