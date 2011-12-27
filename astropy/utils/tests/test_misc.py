# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .. import misc


def test_pkg_finder():
    #note that this also implicitly tests compat.misc._patched_getmodule
    assert misc.find_current_module(0).__name__ == 'astropy.utils.misc'
    assert misc.find_current_module(1).__name__ == 'astropy.utils.tests.test_misc'
    assert misc.find_current_module(0,True).__name__ == 'astropy.utils.tests.test_misc'

class ToBePickled(object):
    def __init__(self,item):
        self.item = item

    def __eq__(self,other):
        if isinstance(other,ToBePickled):
            return self.item == other.item
        else:
            return False


def test_fnpickling(tmpdir):
    fn1 = str(tmpdir.join('test1.pickle'))
    fn2 = str(tmpdir.join('test2.pickle'))
    fn3 = str(tmpdir.join('test3.pickle'))

    obj1 = 'astring'
    misc.fnpickle(obj1, fn1)
    res = misc.fnunpickle(fn1)
    assert obj1 == res

    #try without cPickle and the oldest protocol
    misc.fnpickle(obj1, fn1, usecPickle=False, protocol=0)
    res = misc.fnunpickle(fn1, usecPickle=False)
    assert obj1 == res

    #more complex object
    obj2 = ToBePickled(obj1)
    misc.fnpickle(obj2, fn2)
    res = misc.fnunpickle(fn2)
    assert res == obj2

    #now try multiples
    obj3 = 328.3432
    obj4 = 'blahblahfoo'
    misc.fnpickle(obj3, fn3)
    misc.fnpickle(obj4, fn3, append=True)
    res = misc.fnunpickle(fn3, number=2)
    assert res[0] == obj3
    assert res[1] == obj4
