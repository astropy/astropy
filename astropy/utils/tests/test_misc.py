# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .. import misc


def test_pkg_finder():
    """
    Tests that the `utils.misc.find_current_module` function works. Note that
    this also implicitly tests compat.misc._patched_getmodule
    """
    assert misc.find_current_module(0).__name__ == 'astropy.utils.misc'
    assert misc.find_current_module(1).__name__ == 'astropy.utils.tests.test_misc'
    assert misc.find_current_module(0,True).__name__ == 'astropy.utils.tests.test_misc'

def test_fnpickling_simple(tmpdir):
    """
    Tests the `utils.misc.fnpickle` and `utils.misc.fnupickle` functions' basic
    operation by pickling and unpickling a string, using both a filename and a
    file.
    """
    fn = str(tmpdir.join('test1.pickle'))

    obj1 = 'astring'
    misc.fnpickle(obj1, fn)
    res = misc.fnunpickle(fn)
    assert obj1 == res

    #try without cPickle
    misc.fnpickle(obj1, fn, usecPickle=False)
    res = misc.fnunpickle(fn, usecPickle=False)
    assert obj1 == res

    #now try with a file-like object instead of a string
    with open(fn,'wb') as f:
        misc.fnpickle(obj1, f)
    with open(fn,'rb') as f:
        res = misc.fnunpickle(f)
        assert obj1 == res

    #same without cPickle
    with open(fn,'wb') as f:
        misc.fnpickle(obj1, f, usecPickle=False)
    with open(fn,'rb') as f:
        res = misc.fnunpickle(f, usecPickle=False)
        assert obj1 == res



class ToBePickled(object):
    def __init__(self,item):
        self.item = item

    def __eq__(self,other):
        if isinstance(other,ToBePickled):
            return self.item == other.item
        else:
            return False

def test_fnpickling_class(tmpdir):
    """
    Tests the `utils.misc.fnpickle` and `utils.misc.fnupickle` functions'
    ability to pickle and unpickle custom classes.
    """
    fn = str(tmpdir.join('test2.pickle'))

    obj1 = 'astring'
    obj2 = ToBePickled(obj1)
    misc.fnpickle(obj2, fn)
    res = misc.fnunpickle(fn)
    assert res == obj2

def test_fnpickling_protocol(tmpdir):
    """
    Tests the `utils.misc.fnpickle` and `utils.misc.fnupickle` functions'
    ability to pickle and unpickle pickle files from all protcols.
    """
    import pickle

    obj1 = 'astring'
    obj2 = ToBePickled(obj1)

    for p in range(pickle.HIGHEST_PROTOCOL+1):
        fn = str(tmpdir.join('testp%i.pickle'%p))
        misc.fnpickle(obj2, fn, protocol=p)
        res = misc.fnunpickle(fn)
        assert res == obj2


def test_fnpickling_many(tmpdir):
    """
    Tests the `utils.misc.fnpickle` and `utils.misc.fnupickle` functions'
    ability to pickle and unpickle multiple objects from a single file.
    """
    from pytest import raises

    fn = str(tmpdir.join('test3.pickle'))

    #now try multiples
    obj3 = 328.3432
    obj4 = 'blahblahfoo'
    misc.fnpickle(obj3, fn)
    misc.fnpickle(obj4, fn, append=True)

    res = misc.fnunpickle(fn, number=-1)
    assert len(res) == 2
    assert res[0] == obj3
    assert res[1] == obj4

    misc.fnpickle(obj4, fn, append=True)
    res = misc.fnunpickle(fn, number=2)
    assert len(res) == 2

    with raises(EOFError):
        misc.fnunpickle(fn, number=5)
