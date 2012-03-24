import numpy as np

import timeit

if __name__ == "__main__":
    for ndims in [1,2,3]:
        print "\n%i-dimensional arrays" % ndims
        print " ".join(["%17s" % n for n in ("n","convolve","convolve_fft","convolve_fftw")])

        for ii in xrange(3,15-ndims*3):
            #array = np.random.random([2**ii]*ndims)
            setup=("from astropy.nddata.convolution.convolve import convolve; "+
                "from astropy.nddata.convolution.convolve_fft import convolve_fft; "+
                "from astropy.nddata.convolution.make_kernel import make_kernel; "+
                "import numpy as np; "+
                "array = np.random.random([%i]*%i); "+
                "kernel = make_kernel([%i]*%i,3,force_odd=True)") % (2**ii,ndims,2**ii,ndims)

            #print array.mean(),array.std()
            #print [ffttype for ffttype in ('sp','np','w')]
            #print min(timeit.Timer(stmt="test_ffts.%sfft2(array)" % 'sp',setup=setup).repeat(3,100))

            print "%16i:" % (int(2**ii)) + \
                    "".join(
                        ["%17f" % (min(timeit.Timer(stmt="convolve%s(array,kernel%s)" % (ffttype,extra),setup=setup).repeat(3,10)))
                            for ffttype,extra in zip(("","_fft","_fft"),("",",use_numpy_fft=True",""))]
                    )
