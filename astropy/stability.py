# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This lets code that cares about using only stable APIs declare that
intention.

If code does::

    from astropy import stability

and then subsequently imports a subpackage with an unstable API, a
RuntimeError will be raised.
"""

import astropy

astropy.__require_stable_api__ = True
