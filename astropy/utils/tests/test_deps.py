# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .. import dependencies


def test_opdep():
    #this should always work because numpy is an astropy requirement, but it
    #tests the machinery
    @dependencies.requires_optional_dependencies('numpy')
    def makearray1(n):
        from numpy import random

        return random.randn(n)

    makearray1(10)

    @dependencies.requires_optional_dependencies('numpy, numpy')
    def makearray2(n):
        from numpy import random

        return random.randn(n)

    makearray2(10)

    @dependencies.requires_optional_dependencies('numpy', 'numpy')
    def makearray3(n):
        from numpy import random

        return random.randn(n)

    makearray3(10)


def test_find_deps():
	# note that this will fail if any astropy module isn't importable
	dependencies.find_all_optional_dependencies('astropy')