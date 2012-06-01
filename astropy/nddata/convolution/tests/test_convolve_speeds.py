import numpy as np

import timeit

if __name__ == "__main__":
    for ndims in [1,2,3]:
        print "\n%i-dimensional arrays ('n' is the size of the image AND the kernel)" % ndims
        print " ".join(["%17s" % n for n in ("n","convolve","convolve_fft","convolve_fftw")])

        for ii in xrange(3,15-ndims*3):
            #array = np.random.random([2**ii]*ndims)
            setup=("""
from astropy.nddata.convolution.convolve import convolve; 
from astropy.nddata.convolution.convolve import convolve_fft; 
from astropy.nddata.convolution.make_kernel import make_kernel; 
import numpy as np; 
array = np.random.random([%i]*%i); 
kernel = make_kernel([%i]*%i,3,force_odd=True)""") % (2**ii,ndims,2**ii,ndims)


            print "%16i:" % (int(2**ii)) ,

            for ffttype,extra in zip(("","_fftw","_fftnp","_fftsp"),
                    ("","fft_type='fftw'","fft_type='numpy'","fft_type='scipy'")):
                statement = "convolve%s(array,kernel,boundary='fill',%s)" % (ffttype,extra)
                besttime = min(timeit.Timer(stmt=statement,setup=setup).repeat(3,10))
                print "%17f" % (besttime),

            print

"""
RESULTS on a 2011 Mac Air:
1-dimensional arrays
                n          convolve      convolve_fft     convolve_fftw
               8:          0.000316          0.002252          0.004897
              16:          0.000360          0.002071          0.004604
              32:          0.000488          0.002308          0.005645
              64:          0.000761          0.002554          0.007544
             128:          0.000737          0.002603          0.005438
             256:          0.002322          0.002485          0.006206
             512:          0.007402          0.003402          0.006586
            1024:          0.026499          0.004393          0.007358
            2048:          0.114207          0.007395          0.011030

2-dimensional arrays
                n          convolve      convolve_fft     convolve_fftw
               8:          0.000828          0.003175          0.006379
              16:          0.003191          0.004242          0.006543
              32:          0.044408          0.007972          0.009636
              64:          0.775630          0.027424          0.035715
             128:         12.601176          0.139741          0.138202
             256:        246.371845          0.822725          0.773616

3-dimensional arrays
                n          convolve      convolve_fft     convolve_fftw
               8:          0.012053          0.009994          0.010265
              16:          0.890653          0.071099          0.045073
              32:         61.782118          0.973795          0.654251    
"""
