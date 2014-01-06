# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from ...tests.helper import pytest
from ...extern.six.moves import cPickle
from ... import convolution as conv


def check_pickling_recovery(original, testfile):
    # Try to pickle an object. If successful, make sure 
    # the object's attributes survived pickling and unpickling.
    original_dict = original.__dict__
    with open(testfile, 'w') as f:
        cPickle.dump(original, f)

    with open(testfile, 'r') as f:
        unpickled = cPickle.load(f)
    unpickled_dict = unpickled.__dict__
    
    for original_key in original_dict:
        assert unpickled_dict.has_key(original_key),\
          "Did not pickle {}".format(original_key)
    

@pytest.mark.parametrize("name,args,kwargs", 
                         [(conv.CustomKernel, [], 
                           {'array':np.random.rand(15)}),
                          (conv.Gaussian1DKernel, [1.0],
                           {'x_size':5}),
                          (conv.Gaussian2DKernel, [1.0],
                           {'x_size':5, 'y_size':5}),
                         ])
def test_simple_object(tmpdir, name, args, kwargs):
    # Tests easily instantiated objects
    testfile = str(tmpdir.join('testfile'))

    original = name(*args, **kwargs)
    check_pickling_recovery(original, testfile)
    
