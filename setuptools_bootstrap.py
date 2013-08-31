# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Pre-ez_setup bootstrap module to ensure that either distribute or setuptools >=
0.7 is used (over pre-distribute setuptools) if it is available on the path;
otherwise the latest setuptools will be downloaded and bootstrapped with
``ez_setup.py``.
"""

import sys
import imp

try:
    import pkg_resources
    _setuptools_req = pkg_resources.Requirement.parse('setuptools>=0.7')
    # This may raise a DistributionNotFound in which case no version of
    # setuptools or distribute is properly instlaled
    _setuptools = pkg_resources.get_distribution('setuptools')
    if _setuptools not in _setuptools_req:
        # Older version of setuptools; check if we have distribute; again if
        # this results in DistributionNotFound we want to give up
        _distribute = pkg_resources.get_distribution('distribute')
        if _setuptools != _distribute:
            # It's possible on some pathological systems to have an old version
            # of setuptools and distribute on sys.path simultaneously; make
            # sure distribute is the one that's used
            sys.path.insert(1, _distribute.location)
            _distribute.activate()
            imp.reload(pkg_resources)
except:
    # There are several types of exceptions that can occur here; if all else
    # fails bootstrap and use the bootstrapped version
    from ez_setup import use_setuptools
    use_setuptools()
