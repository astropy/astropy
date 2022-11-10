# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.utils import collections


def test_homogeneous_list():
    l = collections.HomogeneousList(int)
    with pytest.raises(TypeError):
        l.append(5.0)


def test_homogeneous_list2():
    l = collections.HomogeneousList(int)
    with pytest.raises(TypeError):
        l.extend([5.0])


def test_homogeneous_list3():
    l = collections.HomogeneousList(int)
    l.append(5)
    assert l == [5]


def test_homogeneous_list4():
    l = collections.HomogeneousList(int)
    l.extend([5])
    assert l == [5]


def test_homogeneous_list5():
    l = collections.HomogeneousList(int, [1, 2, 3])
    with pytest.raises(TypeError):
        l[1] = 5.0


def test_homogeneous_list_setitem_works():
    l = collections.HomogeneousList(int, [1, 2, 3])
    l[1] = 5
    assert l == [1, 5, 3]


def test_homogeneous_list_setitem_works_with_slice():
    l = collections.HomogeneousList(int, [1, 2, 3])
    l[0:1] = [10, 20, 30]
    assert l == [10, 20, 30, 2, 3]

    l[:] = [5, 4, 3]
    assert l == [5, 4, 3]

    l[::2] = [2, 1]
    assert l == [2, 4, 1]


def test_homogeneous_list_init_got_invalid_type():
    with pytest.raises(TypeError):
        collections.HomogeneousList(int, [1, 2.0, 3])


def test_homogeneous_list_works_with_generators():
    hl = collections.HomogeneousList(int, (i for i in range(3)))
    assert hl == [0, 1, 2]

    hl = collections.HomogeneousList(int)
    hl.extend(i for i in range(3))
    assert hl == [0, 1, 2]

    hl = collections.HomogeneousList(int)
    hl[0:1] = (i for i in range(3))
    assert hl == [0, 1, 2]

    hl = collections.HomogeneousList(int)
    hl += (i for i in range(3))
    assert hl == [0, 1, 2]
