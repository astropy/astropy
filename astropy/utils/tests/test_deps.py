# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .. import dependencies

def test_opdep():
    #this should always work because numpy is an astropy requirement, but it
    #tests the machinery
    @dependencies.requires_optional_dependencies('numpy')
    def makearray(n):
        from numpy import random

        return random.randn(n)

    makearray(10)
