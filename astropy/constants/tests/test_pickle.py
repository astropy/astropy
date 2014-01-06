# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from ...tests.helper import pytest
from ...extern.six.moves import cPickle
from ... import constants as const


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
    

def test_new_constant(tmpdir):
    testfile = str(tmpdir.join('testfile'))

    original = const.Constant('h_fake', "Not Planck", 
                              0.0, 'J s', 0.0, 'fakeref',
                              system='si')
    check_pickling_recovery(original, testfile)
    

def test_existing_constant(tmpdir):
    testfile = str(tmpdir.join('testfile'))

    original = const.h
    check_pickling_recovery(original, testfile)    


def test_existing_em_constant(tmpdir):
    testfile = str(tmpdir.join('testfile'))

    original = const.e
    check_pickling_recovery(original, testfile)  
