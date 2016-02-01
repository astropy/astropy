from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""
The modules in the accuracy testing subpackage are primarily intended for
comparison with "known-good" (or at least "known-familiar") datasets. More
basic functionality and sanity checks are in the main ``coordinates/tests``
testing modules.
"""

N_ACCURACY_TESTS = 10  # the number of samples to use per accuracy test
