# Licensed under a 3-clause BSD style license - see LICENSE.rst

import importlib
import numpy as np

for erfa_wrapper_name in ['cython_numpy', 'cython_numpy_auto']:

    erfa_module_name = '.'.join(['erfa', erfa_wrapper_name])
    print("Testing with {0}...".format(erfa_module_name))

    erfa = importlib.import_module(erfa_module_name)

    jd = np.linspace(2456855.5, 2456855.5+1.0/24.0/60.0, 60*2+1)
    ra  = np.linspace(0.0,np.pi*2.0,5)
    dec = np.linspace(-np.pi/2.0,+np.pi/2.0,4)

    aob, zob, hob, dob, rob, eo = erfa.atco13(0.0,0.0,0.0,0.0,0.0,0.0,jd,0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
    print(aob.shape)

    aob, zob, hob, dob, rob, eo = erfa.atco13(0.0,0.0,0.0,0.0,0.0,0.0,jd[0],0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
    print(aob.shape)

    aob, zob, hob, dob, rob, eo = erfa.atco13(ra[:,None,None],dec[None,:,None],0.0,0.0,0.0,0.0,jd[None,None,:],0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
    print(aob.shape)

    iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd, 0.0)
    print(iy.shape, ihmsf.shape, ihmsf.dtype)
    print(ihmsf)

    iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd[0], 0.0)
    print(iy.shape, ihmsf.shape, ihmsf.dtype)
    print(ihmsf)

    astrom = np.zeros([2],dtype=erfa.dt_eraASTROM)
    theta = np.arange(0,10.0)
    print(theta.shape)
    print(astrom.shape)
    astrom = erfa.aper(theta[:,None], astrom[None,:])
    print(astrom)
