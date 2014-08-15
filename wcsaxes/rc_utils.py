"""
This is a backport of the rc_context class from matplotlib 1.2.
"""

from matplotlib import rcParams

try:
    from matplotlib import rc_context
except ImportError:
    class rc_context(object):
        """
        Return a context manager for managing rc settings.

        This allows one to do::

            with mpl.rc_context(fname='screen.rc'):
                plt.plot(x, a)
                with mpl.rc_context(fname='print.rc'):
                    plt.plot(x, b)
                plt.plot(x, c)

        The 'a' vs 'x' and 'c' vs 'x' plots would have settings from
        'screen.rc', while the 'b' vs 'x' plot would have settings from
        'print.rc'.

        A dictionary can also be passed to the context manager::

            with mpl.rc_context(rc={'text.usetex': True}, fname='screen.rc'):
                plt.plot(x, a)

        The 'rc' dictionary takes precedence over the settings loaded from
        'fname'.  Passing a dictionary only is also valid.
        """

        def __init__(self, rc=None, fname=None):
            self.rcdict = rc
            self.fname = fname
            self._rcparams = rcParams.copy()
            if self.fname:
                rc_file(self.fname)
            if self.rcdict:
                rcParams.update(self.rcdict)

        def __enter__(self):
            return self

        def __exit__(self, type, value, tb):
            rcParams.update(self._rcparams)
