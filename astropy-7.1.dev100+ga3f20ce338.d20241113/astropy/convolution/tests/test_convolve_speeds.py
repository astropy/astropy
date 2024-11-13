# Licensed under a 3-clause BSD style license - see LICENSE.rst

import timeit

# largest image size to use for "linear" and fft convolutions
max_exponents_linear = {1: 15, 2: 7, 3: 5}
max_exponents_fft = {1: 15, 2: 10, 3: 7}

if __name__ == "__main__":
    for ndims in [1, 2, 3]:
        print(
            f"\n{ndims}-dimensional arrays ('n' is the size of the image AND "
            "the kernel)"
        )
        print(" ".join(["%17s" % n for n in ("n", "convolve", "convolve_fft")]))

        for ii in range(3, max_exponents_fft[ndims]):
            # array = np.random.random([2**ii]*ndims)
            # test ODD sizes too
            if ii < max_exponents_fft[ndims]:
                setup = """
import numpy as np
from astropy.convolution.convolve import convolve
from astropy.convolution.convolve import convolve_fft
array = np.random.random([%i]*%i)
kernel = np.random.random([%i]*%i)""" % (
                    2**ii - 1,
                    ndims,
                    2**ii - 1,
                    ndims,
                )

                print("%16i:" % (int(2**ii - 1)), end=" ")

                if ii <= max_exponents_linear[ndims]:
                    for convolve_type, extra in zip(
                        ("", "_fft"), ("", "fft_pad=False")
                    ):
                        statement = (
                            f"convolve{convolve_type}(array, kernel, "
                            f"boundary='fill', {extra})"
                        )
                        besttime = min(
                            timeit.Timer(stmt=statement, setup=setup).repeat(3, 10)
                        )
                        print(f"{besttime:17f}", end=" ")
                else:
                    print("%17s" % "skipped", end=" ")
                    statement = "convolve_fft(array, kernel, boundary='fill')"
                    besttime = min(
                        timeit.Timer(stmt=statement, setup=setup).repeat(3, 10)
                    )
                    print(f"{besttime:17f}", end=" ")

                print()

            setup = """
import numpy as np
from astropy.convolution.convolve import convolve
from astropy.convolution.convolve import convolve_fft
array = np.random.random([%i]*%i)
kernel = np.random.random([%i]*%i)""" % (
                2**ii,
                ndims,
                2**ii,
                ndims,
            )

            print("%16i:" % (int(2**ii)), end=" ")

            if ii <= max_exponents_linear[ndims]:
                for convolve_type in (
                    "",
                    "_fft",
                ):
                    # convolve doesn't allow even-sized kernels
                    if convolve_type == "":
                        print("%17s" % "-", end=" ")
                    else:
                        statement = (
                            f"convolve{convolve_type}(array, kernel, boundary='fill')"
                        )
                        besttime = min(
                            timeit.Timer(stmt=statement, setup=setup).repeat(3, 10)
                        )
                        print(f"{besttime:17f}", end=" ")
            else:
                print("%17s" % "skipped", end=" ")
                statement = "convolve_fft(array, kernel, boundary='fill')"
                besttime = min(timeit.Timer(stmt=statement, setup=setup).repeat(3, 10))
                print(f"{besttime:17f}", end=" ")

            print()

"""
Unfortunately, these tests are pretty strongly inconclusive
NOTE: Runtime has units seconds and represents wall clock time.

RESULTS on a late 2013 Mac Pro:
3.5 GHz 6-Core Intel Xeon E5
32 GB 1866 MHz DDR3 ECC
Python 3.5.4 :: Anaconda custom (x86_64)
clang version 6.0.0 (tags/RELEASE_600/final)
llvm-opnemp r327556 | grokos | 2018-03-14 15:11:36 -0400 (Wed, 14 Mar 2018)

With OpenMP (hyperthreaded 12procs), convolve() only:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.002895          0.007321
              15:          0.002684          0.008028
              31:          0.002733          0.008684
              63:          0.002728          0.009127
             127:          0.002851          0.012659
             255:          0.002835          0.010550
             511:          0.003051          0.017137
            1023:          0.004042          0.019384
            2047:          0.007371          0.049246
            4095:          0.021903          0.039821
            8191:          0.067098          8.335749
           16383:          0.256072          0.272165

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.002696          0.014745
              15:          0.002839          0.014826
              31:          0.004286          0.045167
              63:          0.022941          0.063715
             127:          0.325557          0.925577
             255:           skipped          0.694621
             511:           skipped          3.734946

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.003502          0.033121
               8:          0.003407          0.030351
              15:          0.026338          0.062235
              31:          1.239503          1.586930
              63:           skipped         10.792675

With OpenMP but single threaded (n_threads = 1), convolve() only:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.001754          0.004687
              15:          0.001706          0.005133
              31:          0.001744          0.005381
              63:          0.001725          0.005582
             127:          0.001801          0.007405
             255:          0.002262          0.006528
             511:          0.003866          0.009913
            1023:          0.009820          0.011511
            2047:          0.034707          0.028171
            4095:          0.132908          0.024133
            8191:          0.527692          8.311933
           16383:          2.103046          0.269368

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.001734          0.009458
              15:          0.002336          0.010310
              31:          0.009123          0.025427
              63:          0.126701          0.040610
             127:          2.126114          0.926549
             255:           skipped          0.690896
             511:           skipped          3.756475

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.002822          0.019498
              15:          0.096008          0.063744
              31:          7.373533          1.578913
              63:           skipped         10.811530

RESULTS on a 2011 Mac Air:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000408          0.002334          0.005571          0.002677
              15:          0.000361          0.002491          0.005648          0.002678
              31:          0.000535          0.002450          0.005988          0.002880
              63:          0.000509          0.002876          0.008003          0.002981
             127:          0.000801          0.004080          0.008513          0.003932
             255:          0.002453          0.003111          0.007518          0.003564
             511:          0.008394          0.006224          0.010247          0.005991
            1023:          0.028741          0.007538          0.009591          0.007696
            2047:          0.106323          0.021575          0.022041          0.020682
            4095:          0.411936          0.021675          0.019761          0.020939
            8191:          1.664517          8.278320          0.073001          7.803563
           16383:          6.654678          0.251661          0.202271          0.222171

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000552          0.003524          0.006667          0.004318
              15:          0.002986          0.005093          0.012941          0.005951
              31:          0.074360          0.033973          0.031800          0.036937
              63:          0.848471          0.057407          0.052192          0.053213
             127:         14.656414          1.005329          0.402113          0.955279
             255:           skipped          1.715546          1.566876          1.745338
             511:           skipped          4.066155          4.303350          3.930661

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.009239          0.012957          0.011957          0.015997
              15:          0.772434          0.075621          0.056711          0.079508
              31:         62.824051          2.295193          1.189505          2.351136
              63:           skipped         11.250225         10.982726         10.585744



On a 2009 Mac Pro:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000360          0.002269          0.004986          0.002476
              15:          0.000364          0.002255          0.005244          0.002471
              31:          0.000385          0.002380          0.005422          0.002588
              63:          0.000474          0.002407          0.005392          0.002637
             127:          0.000752          0.004122          0.007827          0.003966
             255:          0.004316          0.003258          0.006566          0.003324
             511:          0.011517          0.007158          0.009898          0.006238
            1023:          0.034105          0.009211          0.009468          0.008260
            2047:          0.113620          0.028097          0.020662          0.021603
            4095:          0.403373          0.023211          0.018767          0.020065
            8191:          1.519329          8.454573          0.211436          7.212381
           16383:          5.887481          0.317428          0.153344          0.237119

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000474          0.003470          0.006131          0.003503
              15:          0.002011          0.004481          0.007825          0.004496
              31:          0.027291          0.019433          0.014841          0.018034
              63:          0.445680          0.038171          0.026753          0.037404
             127:          7.003774          0.925921          0.282591          0.762671
             255:           skipped          0.804682          0.708849          0.869368
             511:           skipped          3.643626          3.687562          4.584770

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.004520          0.011519          0.009464          0.012335
              15:          0.329566          0.060978          0.045495          0.073692
              31:         24.935228          1.654920          0.710509          1.773879
              63:           skipped          8.982771         12.407683         16.900078
"""
