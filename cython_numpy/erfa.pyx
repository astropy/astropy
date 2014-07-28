import numpy as np
cimport numpy as np

np.import_array()

__all__ = ['atco13', 'd2dtf', 'aper']


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
    void eraAper(double theta, eraASTROM *astrom)

cdef struct eraASTROM:
    double pmt
    double eb[3]
    double eh[3]
    double em
    double v[3]
    double bm1 
    double bpn[3][3]
    double along
    double phi
    double xpl
    double ypl
    double sphi
    double cphi
    double diurab   
    double eral     
    double refa
    double refb

dt_eraASTROM = np.dtype([('pmt','d'),
                         ('eb','d',(3,)),
                         ('eh','d',(3,)),
                         ('em','d'),
                         ('v','d',(3,)),
                         ('bm1 ','d'),
                         ('bpn','d',(3,3)),
                         ('along','d'),
                         ('phi','d'),
                         ('xpl','d'),
                         ('ypl','d'),
                         ('sphi','d'),
                         ('cphi','d'),
                         ('diurab','d'),
                         ('eral','d'),
                         ('refa','d'),
                         ('refb','d')], align=True)


#Note: the pattern used here follows https://github.com/cython/cython/wiki/tutorials-numpy#dimensionally-simple-functions



def atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl):
    
    shape = np.broadcast(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl).shape    
    aob_out = np.empty(shape, dtype=np.double)
    zob_out = np.empty(shape, dtype=np.double)
    hob_out = np.empty(shape, dtype=np.double)
    dob_out = np.empty(shape, dtype=np.double)
    rob_out = np.empty(shape, dtype=np.double)
    eo_out  = np.empty(shape, dtype=np.double)
    
    cdef np.broadcast it = np.broadcast(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl, aob_out, zob_out, hob_out, dob_out, rob_out, eo_out)
    
    cdef double _rc
    cdef double _dc
    cdef double _pr
    cdef double _pd
    cdef double _px
    cdef double _rv
    cdef double _utc1
    cdef double _utc2
    cdef double _dut1
    cdef double _elong
    cdef double _phi
    cdef double _hm
    cdef double _xp
    cdef double _yp
    cdef double _phpa
    cdef double _tk
    cdef double _rh
    cdef double _wl
    cdef double *_aob
    cdef double *_zob
    cdef double *_hob
    cdef double *_dob
    cdef double *_rob
    cdef double *_eo
    
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
        _aob   = (<double*>np.PyArray_MultiIter_DATA(it, 18))
        _zob   = (<double*>np.PyArray_MultiIter_DATA(it, 19)) 
        _hob   = (<double*>np.PyArray_MultiIter_DATA(it, 20))
        _dob   = (<double*>np.PyArray_MultiIter_DATA(it, 21))
        _rob   = (<double*>np.PyArray_MultiIter_DATA(it, 22))
        _eo    = (<double*>np.PyArray_MultiIter_DATA(it, 23))
        
        ret = eraAtco13(_rc, _dc, _pr, _pd, _px, _rv, _utc1, _utc2, _dut1, _elong, _phi, _hm, _xp, _yp, _phpa, _tk, _rh, _wl, _aob, _zob, _hob, _dob, _rob, _eo)        
        
        np.PyArray_MultiIter_NEXT(it)
    
    return aob_out, zob_out, hob_out, dob_out, rob_out, eo_out



def d2dtf(scale, ndp, d1, d2):
    
    shape = np.broadcast(scale, ndp, d1, d2).shape
    iy_out    = np.empty(shape, dtype=np.int)
    im_out    = np.empty(shape, dtype=np.int)
    id_out    = np.empty(shape, dtype=np.int)
    ihmsf_out = np.empty(shape, dtype=[('h','i'),('m','i'),('s','i'),('f','i')]) 
    
    cdef np.broadcast it = np.broadcast(scale, ndp, d1, d2, iy_out, im_out, id_out, ihmsf_out)
    
    cdef char *_scale
    cdef int _ndp
    cdef double _d1
    cdef double _d2
    cdef int *_iy
    cdef int *_im
    cdef int *_id
    cdef int *_ihmsf
    
    while np.PyArray_MultiIter_NOTDONE(it):
        
        _scale = (  <char*>np.PyArray_MultiIter_DATA(it,  0))
        _ndp   = (   <int*>np.PyArray_MultiIter_DATA(it,  1))[0]
        _d1    = (<double*>np.PyArray_MultiIter_DATA(it,  2))[0]
        _d2    = (<double*>np.PyArray_MultiIter_DATA(it,  3))[0]
        _iy    = (   <int*>np.PyArray_MultiIter_DATA(it,  4))
        _im    = (   <int*>np.PyArray_MultiIter_DATA(it,  5))
        _id    = (   <int*>np.PyArray_MultiIter_DATA(it,  6))
        _ihmsf = (   <int*>np.PyArray_MultiIter_DATA(it,  7))
        
        ret = eraD2dtf(_scale, _ndp, _d1, _d2, _iy, _im, _id, _ihmsf)
        
        np.PyArray_MultiIter_NEXT(it)
    
    return iy_out, im_out, id_out, ihmsf_out



def aper(theta, astrom):
    
    shape = np.broadcast(theta, astrom).shape
    astrom_out = np.empty(shape, dtype=dt_eraASTROM)
    np.copyto(astrom_out, astrom)
    
    cdef np.broadcast it = np.broadcast(theta, astrom_out)
    
    cdef double _theta
    cdef eraASTROM *_astrom
    
    while np.PyArray_MultiIter_NOTDONE(it):
        
        _theta  = (   <double *>np.PyArray_MultiIter_DATA(it,  0))[0]
        _astrom = (<eraASTROM *>np.PyArray_MultiIter_DATA(it,  1))
        
        eraAper(_theta, _astrom)
        
        np.PyArray_MultiIter_NEXT(it)
    
    return astrom_out