import numpy as np
cimport numpy as np

np.import_array()

__all__ = ['atco13', 'd2dtf']

cdef extern from "erfa.h":
    int eraAtco13(double rc, double dc,
              double pr, double pd, double px, double rv,
              double utc1, double utc2, double dut1,
              double elong, double phi, double hm, double xp, double yp,
              double phpa, double tk, double rh, double wl,
              double *aob, double *zob, double *hob,
              double *dob, double *rob, double *eo)
    int eraD2dtf(const char *scale, int ndp, double d1, double d2,
                 int *iy, int *im, int *id, int ihmsf[4])

#Note: the pattern used here follows https://github.com/cython/cython/wiki/tutorials-numpy#dimensionally-simple-functions


def atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl):
    
    shape = np.broadcast(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl).shape    
    aob = np.empty(shape, dtype=np.double)
    zob = np.empty(shape, dtype=np.double)
    hob = np.empty(shape, dtype=np.double)
    dob = np.empty(shape, dtype=np.double)
    rob = np.empty(shape, dtype=np.double)
    eo  = np.empty(shape, dtype=np.double)
    
    cdef np.broadcast it = np.broadcast(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl, aob, zob, hob, dob, rob, eo)
    
    cdef double _aob
    cdef double _zob
    cdef double _hob
    cdef double _dob
    cdef double _rob
    cdef double _eo
    
    while np.PyArray_MultiIter_NOTDONE(it):
        
        _rc    = (<double*>np.PyArray_MultiIter_DATA(it,  0))[0]
        _dc    = (<double*>np.PyArray_MultiIter_DATA(it,  1))[0]
        _pr    = (<double*>np.PyArray_MultiIter_DATA(it,  2))[0]
        _pd    = (<double*>np.PyArray_MultiIter_DATA(it,  3))[0]
        _px    = (<double*>np.PyArray_MultiIter_DATA(it,  4))[0]
        _rv    = (<double*>np.PyArray_MultiIter_DATA(it,  5))[0]
        _utc1  = (<double*>np.PyArray_MultiIter_DATA(it,  6))[0]
        _utc2  = (<double*>np.PyArray_MultiIter_DATA(it,  7))[0]
        _dut1  = (<double*>np.PyArray_MultiIter_DATA(it,  8))[0]
        _elong = (<double*>np.PyArray_MultiIter_DATA(it,  9))[0]
        _phi   = (<double*>np.PyArray_MultiIter_DATA(it, 10))[0]
        _hm    = (<double*>np.PyArray_MultiIter_DATA(it, 11))[0]
        _xp    = (<double*>np.PyArray_MultiIter_DATA(it, 12))[0]
        _yp    = (<double*>np.PyArray_MultiIter_DATA(it, 13))[0]
        _phpa  = (<double*>np.PyArray_MultiIter_DATA(it, 14))[0]
        _tk    = (<double*>np.PyArray_MultiIter_DATA(it, 15))[0]
        _rh    = (<double*>np.PyArray_MultiIter_DATA(it, 16))[0]
        _wl    = (<double*>np.PyArray_MultiIter_DATA(it, 17))[0]
        
        ret = eraAtco13(_rc, _dc, _pr, _pd, _px, _rv, _utc1, _utc2, _dut1, _elong, _phi, _hm, _xp, _yp, _phpa, _tk, _rh, _wl, &_aob, &_zob, &_hob, &_dob, &_rob, &_eo)
        
        (<double*>np.PyArray_MultiIter_DATA(it, 18))[0] = _aob
        (<double*>np.PyArray_MultiIter_DATA(it, 19))[0] = _zob
        (<double*>np.PyArray_MultiIter_DATA(it, 20))[0] = _hob
        (<double*>np.PyArray_MultiIter_DATA(it, 21))[0] = _dob
        (<double*>np.PyArray_MultiIter_DATA(it, 22))[0] = _rob
        (<double*>np.PyArray_MultiIter_DATA(it, 23))[0] = _eo
        
        np.PyArray_MultiIter_NEXT(it)
    
    return aob, zob, hob, dob, rob, eo

def d2dtf(scale, ndp, d1, d2):
    
    shape = np.broadcast(d1, d2).shape
    iy    = np.empty(shape, dtype=np.int)
    im    = np.empty(shape, dtype=np.int)
    id    = np.empty(shape, dtype=np.int)
    ihmsf = np.empty(shape, dtype=[('h','i'),('m','i'),('s','i'),('f','i')]) 
    
    cdef np.broadcast it = np.broadcast(d1, d2, iy, im, id, ihmsf)
    
    cdef int _iy
    cdef int _im
    cdef int _id
    cdef int _ihmsf[4]
    
    while np.PyArray_MultiIter_NOTDONE(it):
        
        _d1    = (<double*>np.PyArray_MultiIter_DATA(it,  0))[0]
        _d2    = (<double*>np.PyArray_MultiIter_DATA(it,  1))[0]
        
        ret = eraD2dtf(scale, ndp, _d1, _d2, &_iy, &_im, &_id, _ihmsf)
        
        (<int*>np.PyArray_MultiIter_DATA(it, 2))[0] = _iy
        (<int*>np.PyArray_MultiIter_DATA(it, 3))[0] = _im
        (<int*>np.PyArray_MultiIter_DATA(it, 4))[0] = _id
        (<int*>np.PyArray_MultiIter_DATA(it, 5))[0:4] = _ihmsf
        
        np.PyArray_MultiIter_NEXT(it)
    
    return iy, im, id, ihmsf
