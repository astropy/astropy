# Licensed under a 3-clause BSD style license - see LICENSE.rst

def test_nddata_basic():
    from numpy.random import randn,rand
    from numpy import ones,dtype
    from ..nddata import NDData
    from pytest import raises

    with raises(TypeError):
        NDData()  # empty initializer should fail

    nd1 = NDData(randn(10,10))
    assert nd1.shape == (10,10)
    assert nd1.size == 100
    assert nd1.dtype == dtype(float)

    # now one with error and masks
    nd2 = NDData(randn(10,10),rand(10,10)/5,ones((10,10),dtype=bool))

    with raises(ValueError):
        #mismatched error shape
        nd3 = NDData(randn(10,10),rand(10,9))

    with raises(ValueError):
        #mismatched mask shape
        nd4 = NDData(randn(10,10),mask=ones((10,9),dtype=float))


def test_nddata_copy_conv():
    from numpy.random import randn
    from numpy import dtype
    from ..nddata import NDData

    a = randn(10,10)

    #check that copy works as expected
    ndcopy = NDData(a,copy=True)
    ndref = NDData(a,copy=False)
    a[0,0] = 0
    assert ndref.data[0,0] == 0
    assert ndcopy.data[0,0] != 0

    #conversion
    nd3 = NDData([[1,2,3],[4,5,6]])
    assert nd3.size == 6
    assert nd3.dtype == dtype(int)

def test_boolmask():
    from numpy.random import rand,randn
    from numpy import array,arange,all
    from ..nddata import NDData
    import random, string

    a = randn(10,10)
    bm = randn(10,10)>0  # random mask that boolmask should look like

    nd1 = NDData(randn(10,10),mask=~bm) # mask False where valid
    assert all(nd1.boolmask == bm)

    #scalar mask 0 where valid
    im = arange(100).reshape((10,10)) + 1
    im[bm] = 0
    nd2 = NDData(randn(10,10),mask=im)
    assert all(nd2.boolmask == bm)


    #generate a list of random strings
    strlst = []
    for N in rand(100)*9+1:
        sl = [random.choice(string.ascii_letters) for x in range(int(N))]
        strlst.append(''.join(sl))

    #string mask '' where valid
    strmsk = array(strlst).reshape((10,10))
    strmsk[bm] = ''
    nd3 = NDData(randn(10,10),mask=strmsk)
    assert all(nd3.boolmask == bm)
    #do it with unicode too for goo measure
    umsk = array([unicode(s) for s in strlst]).reshape((10,10))
    umsk[bm] = ''
    nd3 = NDData(randn(10,10),mask=umsk)
    assert all(nd3.boolmask == bm)
