# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import raises

from .. import collections


@raises(TypeError)
def test_homogeneous_list():
    l = collections.HomogeneousList(int)
    l.append(5.0)


@raises(TypeError)
def test_homogeneous_list2():
    l = collections.HomogeneousList(int)
    l.extend([5.0])


def test_homogeneous_list3():
    l = collections.HomogeneousList(int)
    l.append(5)


def test_homogeneous_list4():
    l = collections.HomogeneousList(int)
    l.extend([5])


@raises(TypeError)
def test_homogeneous_list5():
    l = collections.HomogeneousList(int, [1, 2, 3])
    l[1] = 5.0
