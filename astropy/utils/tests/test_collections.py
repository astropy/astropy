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


class OwnedList(collections._OwnedMixin, list):
    pass


class ListOwner:
    def __init__(self, data: list):
        self._data = OwnedList(
            data,
            owner_class=self.__class__,
            public_attr_name="data",
            replacements={
                "append": "add_item",
                "extend": "add_item",
            },
        )

    @property
    def data(self):
        return self._data

    def add_item(self, item):
        self._data.append(item, _owned=True)


def test_ListOwner_add_item():
    lo = ListOwner([1, 2])
    lo.add_item(1)


@pytest.mark.parametrize(
    "method, args",
    [
        ("append", (1,)),
        ("clear", ()),
        ("extend", ([3, 4],)),
        ("insert", (1, 999)),
        ("pop", ()),
        ("remove", (1,)),
        ("reverse", ()),
        ("sort", ()),
        ("__iadd__", ([3, 4],)),
        ("__imul__", (2,)),
        ("__setitem__", (0, 999)),
        ("__delitem__", (1,)),
    ],
)
def test_direct_mutation_owned_list(method, args):
    lo = ListOwner([1, 2])
    fun = getattr(lo.data, method)
    with pytest.warns(
        DeprecationWarning,
        match="Direct mutations of ListOwner.data",
    ):
        fun(*args)
