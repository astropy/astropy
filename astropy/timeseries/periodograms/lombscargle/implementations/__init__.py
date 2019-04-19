"""Various implementations of the Lomb-Scargle Periodogram"""

from .main import lombscargle, available_methods
from .chi2_impl import lombscargle_chi2
from .scipy_impl import lombscargle_scipy
from .slow_impl import lombscargle_slow
from .fast_impl import lombscargle_fast
from .fastchi2_impl import lombscargle_fastchi2
