# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.io.misc import fnpickle, fnunpickle


def test_fnpickling_simple(tmpdir):
    """
    Tests the `fnpickle` and `fnupickle` functions' basic operation by
    pickling and unpickling a string, using both a filename and a
    file.
    """
    fn = str(tmpdir.join('test1.pickle'))

    obj1 = 'astring'
    fnpickle(obj1, fn)
    res = fnunpickle(fn, 0)
    assert obj1 == res

    # now try with a file-like object instead of a string
    with open(fn, 'wb') as f:
        fnpickle(obj1, f)
    with open(fn, 'rb') as f:
        res = fnunpickle(f)
        assert obj1 == res


class ToBePickled:
    def __init__(self, item):
        self.item = item

    def __eq__(self, other):
        if isinstance(other, ToBePickled):
            return self.item == other.item
        else:
            return False


def test_fnpickling_class(tmpdir):
    """
    Tests the `fnpickle` and `fnupickle` functions' ability to pickle
    and unpickle custom classes.
    """
    fn = str(tmpdir.join('test2.pickle'))

    obj1 = 'astring'
    obj2 = ToBePickled(obj1)
    fnpickle(obj2, fn)
    res = fnunpickle(fn)
    assert res == obj2


def test_fnpickling_protocol(tmpdir):
    """
    Tests the `fnpickle` and `fnupickle` functions' ability to pickle
    and unpickle pickle files from all protcols.
    """
    import pickle

    obj1 = 'astring'
    obj2 = ToBePickled(obj1)

    for p in range(pickle.HIGHEST_PROTOCOL + 1):
        fn = str(tmpdir.join(f'testp{p}.pickle'))
        fnpickle(obj2, fn, protocol=p)
        res = fnunpickle(fn)
        assert res == obj2


def test_fnpickling_many(tmpdir):
    """
    Tests the `fnpickle` and `fnupickle` functions' ability to pickle
    and unpickle multiple objects from a single file.
    """

    fn = str(tmpdir.join('test3.pickle'))

    # now try multiples
    obj3 = 328.3432
    obj4 = 'blahblahfoo'
    fnpickle(obj3, fn)
    fnpickle(obj4, fn, append=True)

    res = fnunpickle(fn, number=-1)
    assert len(res) == 2
    assert res[0] == obj3
    assert res[1] == obj4

    fnpickle(obj4, fn, append=True)
    res = fnunpickle(fn, number=2)
    assert len(res) == 2

    with pytest.raises(EOFError):
        fnunpickle(fn, number=5)
