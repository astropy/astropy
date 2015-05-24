# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# TEST_UNICODE_LITERALS

import numpy as np

from ...extern import six
from ... import units as u
from ... import time

def test_info(table_types):
    """
    Test the info() method of printing a table summary
    """
    a = np.array([1, 2, 3], dtype='int32')
    b = np.array([1, 2, 3], dtype='float32')
    c = np.array(['a', 'c', 'e'], dtype='|S1')
    t = table_types.Table([a, b, c], names=['a', 'b', 'c'])
    masked = 'masked=True ' if t.masked else ''
    string8 = 'string8' if six.PY2 else ' bytes8'

    # Minimal output for a typical table
    out = six.moves.cStringIO()
    t.info(out=out)
    assert out.getvalue().splitlines() == ['<{0} {1}length=3>'.format(t.__class__.__name__, masked),
                                           'name  dtype ',
                                           '---- -------',
                                           '   a   int32',
                                           '   b float32',
                                           '   c {0}'.format(string8)]

    # All output fields including a mixin column
    t['d'] = [1,2,3] * u.m
    t['d'].description = 'description'
    t['a'].format = '%02d'
    t['e'] = time.Time([1,2,3], format='cxcsec')
    out = six.moves.cStringIO()
    t.info(out=out)
    assert out.getvalue().splitlines() == ['<{0} {1}length=3>'.format(t.__class__.__name__, masked),
                                           'name  dtype  unit format description class',
                                           '---- ------- ---- ------ ----------- -----',
                                           '   a   int32        %02d                  ',
                                           '   b float32                              ',
                                           '   c {0}                              '.format(string8),
                                           '   d float64    m        description      ',
                                           '   e  object                          Time']
