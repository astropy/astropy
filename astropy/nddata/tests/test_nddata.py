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

    with raises(TypeError):
        #mask not bool
        nd4 = NDData(randn(10,10),mask=ones((10,10),dtype=float))




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
