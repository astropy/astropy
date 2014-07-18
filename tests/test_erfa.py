import numpy as np
import erfa.cython_numpy as erfa

jd = np.linspace(2456855.5, 2456855.5+1.0/24.0, 3600*2+1)

aob, zob, hob, dob, rob, eo = erfa.atco13(0.0,0.0,0.0,0.0,0.0,0.0,jd,0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
print(aob.shape)

aob, zob, hob, dob, rob, eo = erfa.atco13(0.0,0.0,0.0,0.0,0.0,0.0,jd[0],0.0,0.0,0.0,np.pi/4.0,0.0,0.0,0.0,1014.0,0.0,0.0,0.5)
print(aob.shape)

iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd, 0.0)
print(iy.shape, ihmsf.shape, ihmsf.dtype)
print(ihmsf)

iy, im, id, ihmsf = erfa.d2dtf("UTC", 3, jd[0], 0.0)
print(iy.shape, ihmsf.shape, ihmsf.dtype)
print(ihmsf)

