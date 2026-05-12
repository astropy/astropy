from contextlib import ExitStack, contextmanager

from astropy.visualization.time import time_support
from astropy.visualization.units import quantity_support

__all__ = ["astropy_types_support"]


@contextmanager
def astropy_types_support(*, quantity_support_kwargs=None, time_support_kwargs=None):
    """
    Enable support for plotting `astropy.units.Quantity` and `astropy.time.Time` instances in
    matplotlib.

    It can be used as a decorator or with a ``with`` statement.

    Examples
    --------

    .. plot::
        :include-source:

        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
        import astropy.units as u
        from astropy.time import Time
        from astropy.visualization.astropy_types_support import astropy_types_support

        @astropy_types_support()
        def plot_example():
            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot()
            times = Time(["2000-01-01", "2000-01-02", "2000-01-03"])
            ax.plot([1, 2, 3] * u.m, times)
            canvas.draw()
            return fig

        with astropy_types_support():  # doctest: +IGNORE_OUTPUT
            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot()
            times = Time(["2000-01-01", "2000-01-02", "2000-01-03"])
            ax.plot([1, 2, 3] * u.m, times)
            canvas.draw()

    """
    with ExitStack() as stack:
        stack.enter_context(quantity_support(**(quantity_support_kwargs or {})))
        stack.enter_context(time_support(**(time_support_kwargs or {})))
        yield
