# Licensed under a 3-clause BSD style license - see LICENSE.rst

import timeit

import numpy as np  # pylint: disable=W0611


# largest image size to use for "linear" and fft convolutions
max_exponents_linear = {1: 15, 2: 7, 3: 5}
max_exponents_fft = {1: 15, 2: 10, 3: 7}

if __name__ == "__main__":
    for ndims in [1, 2, 3]:
        print("\n{}-dimensional arrays ('n' is the size of the image AND "
              "the kernel)".format(ndims))
        print(" ".join(["%17s" % n for n in ("n", "convolve", "convolve_fft")]))

        for ii in range(3, max_exponents_fft[ndims]):
            # array = np.random.random([2**ii]*ndims)
            # test ODD sizes too
            if ii < max_exponents_fft[ndims]:
                setup = ("""
import numpy as np
from astropy.convolution.convolve import convolve
from astropy.convolution.convolve import convolve_fft
array = np.random.random([%i]*%i)
kernel = np.random.random([%i]*%i)""") % (2 ** ii - 1, ndims, 2 ** ii - 1, ndims)

                print("%16i:" % (int(2 ** ii - 1)), end=' ')

                if ii <= max_exponents_linear[ndims]:
                    for convolve_type, extra in zip(("", "_fft"),
                                              ("", "fft_pad=False")):
                        statement = "convolve{}(array, kernel, boundary='fill', {})".format(convolve_type, extra)
                        besttime = min(timeit.Timer(stmt=statement, setup=setup).repeat(3, 10))
                        print("%17f" % (besttime), end=' ')
                else:
                    print("%17s" % "skipped", end=' ')
                    statement = "convolve_fft(array, kernel, boundary='fill')"
                    besttime = min(timeit.Timer(stmt=statement, setup=setup).repeat(3, 10))
                    print("%17f" % (besttime), end=' ')

                print()

            setup = ("""
import numpy as np
from astropy.convolution.convolve import convolve
from astropy.convolution.convolve import convolve_fft
array = np.random.random([%i]*%i)
kernel = np.random.random([%i]*%i)""") % (2 ** ii - 1, ndims, 2 ** ii - 1, ndims)

            print("%16i:" % (int(2 ** ii)), end=' ')

            if ii <= max_exponents_linear[ndims]:
                for convolve_type in ("", "_fft"):
                    statement = "convolve{}(array, kernel, boundary='fill')".format(convolve_type)
                    besttime = min(timeit.Timer(stmt=statement, setup=setup).repeat(3, 10))
                    print("%17f" % (besttime), end=' ')
            else:
                print("%17s" % "skipped", end=' ')
                statement = "convolve_fft(array, kernel, boundary='fill')"
                besttime = min(timeit.Timer(stmt=statement, setup=setup).repeat(3, 10))
                print("%17f" % (besttime), end=' ')

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
               8:          0.002737          0.008204
              15:          0.002684          0.008028
              16:          0.002680          0.008682
              31:          0.002733          0.008684
              32:          0.002927          0.009021
              63:          0.002728          0.009127
              64:          0.003757          0.009604
             127:          0.002851          0.012659
             128:          0.002784          0.009814
             255:          0.002835          0.010550
             256:          0.002886          0.010719
             511:          0.003051          0.017137
             512:          0.003176          0.011805
            1023:          0.004042          0.019384
            1024:          0.004074          0.014280
            2047:          0.007371          0.049246
            2048:          0.007535          0.019528
            4095:          0.021903          0.039821
            4096:          0.022138          0.031596
            8191:          0.067098          8.335749
            8192:          0.067217          0.055224
           16383:          0.256072          0.272165
           16384:          0.257656          0.063409

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.002696          0.014745
               8:          0.002796          0.014153
              15:          0.002839          0.014826
              16:          0.002933          0.017755
              31:          0.004286          0.045167
              32:          0.005059          0.024783
              63:          0.022941          0.063715
              64:          0.022747          0.056179
             127:          0.325557          0.925577
             128:          0.325050          0.123447
             255:           skipped          0.694621
             256:           skipped          0.688915
             511:           skipped          3.734946
             512:           skipped          3.735681

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.003502          0.033121
               8:          0.003407          0.030351
              15:          0.026338          0.062235
              16:          0.026138          0.071622
              31:          1.239503          1.586930
              32:          1.237914          0.728493
              63:           skipped         10.792675
              64:           skipped         10.772493

With OpenMP but single threaded (n_threads = 1), convolve() only:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.001754          0.004687
               8:          0.001727          0.005307
              15:          0.001706          0.005133
              16:          0.001671          0.005348
              31:          0.001744          0.005381
              32:          0.001674          0.005603
              63:          0.001725          0.005582
              64:          0.001709          0.005800
             127:          0.001801          0.007405
             128:          0.001813          0.006011
             255:          0.002262          0.006528
             256:          0.002214          0.006413
             511:          0.003866          0.009913
             512:          0.003712          0.007204
            1023:          0.009820          0.011511
            1024:          0.009815          0.008555
            2047:          0.034707          0.028171
            2048:          0.034690          0.011243
            4095:          0.132908          0.024133
            4096:          0.133138          0.018282
            8191:          0.527692          8.311933
            8192:          0.531433          0.031895
           16383:          2.103046          0.269368
           16384:          2.115349          0.063068

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.001734          0.009458
               8:          0.001749          0.009504
              15:          0.002336          0.010310
              16:          0.002156          0.010867
              31:          0.009123          0.025427
              32:          0.009153          0.014798
              63:          0.126701          0.040610
              64:          0.126231          0.033055
             127:          2.126114          0.926549
             128:          2.129131          0.121623
             255:           skipped          0.690896
             256:           skipped          0.689538
             511:           skipped          3.756475
             512:           skipped          3.742207

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve      convolve_fft
               7:          0.002822          0.019498
               8:          0.002784          0.018316
              15:          0.096008          0.063744
              16:          0.095308          0.071917
              31:          7.373533          1.578913
              32:          7.345962          0.734147
              63:           skipped         10.811530
              64:           skipped         10.807303

RESULTS on a 2011 Mac Air:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000408          0.002334          0.005571          0.002677
               8:          0.000399          0.002818          0.006505          0.003094
              15:          0.000361          0.002491          0.005648          0.002678
              16:          0.000371          0.002997          0.005983          0.003036
              31:          0.000535          0.002450          0.005988          0.002880
              32:          0.000452          0.002618          0.007102          0.004366
              63:          0.000509          0.002876          0.008003          0.002981
              64:          0.000453          0.002706          0.005520          0.003049
             127:          0.000801          0.004080          0.008513          0.003932
             128:          0.000749          0.003332          0.006236          0.003159
             255:          0.002453          0.003111          0.007518          0.003564
             256:          0.002478          0.003341          0.006325          0.004290
             511:          0.008394          0.006224          0.010247          0.005991
             512:          0.007934          0.003764          0.006840          0.004106
            1023:          0.028741          0.007538          0.009591          0.007696
            1024:          0.027900          0.004871          0.009628          0.005118
            2047:          0.106323          0.021575          0.022041          0.020682
            2048:          0.108916          0.008107          0.011049          0.007596
            4095:          0.411936          0.021675          0.019761          0.020939
            4096:          0.408992          0.018870          0.016663          0.012890
            8191:          1.664517          8.278320          0.073001          7.803563
            8192:          1.657573          0.037967          0.034227          0.028390
           16383:          6.654678          0.251661          0.202271          0.222171
           16384:          6.611977          0.073630          0.067616          0.055591

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000552          0.003524          0.006667          0.004318
               8:          0.000646          0.004443          0.007354          0.003958
              15:          0.002986          0.005093          0.012941          0.005951
              16:          0.003549          0.005688          0.008818          0.006300
              31:          0.074360          0.033973          0.031800          0.036937
              32:          0.077338          0.017708          0.025637          0.011883
              63:          0.848471          0.057407          0.052192          0.053213
              64:          0.773061          0.029657          0.033409          0.028230
             127:         14.656414          1.005329          0.402113          0.955279
             128:         15.867796          0.266233          0.268551          0.237930
             255:           skipped          1.715546          1.566876          1.745338
             256:           skipped          1.515616          1.268220          1.036881
             511:           skipped          4.066155          4.303350          3.930661
             512:           skipped          3.976139          4.337525          3.968935

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.009239          0.012957          0.011957          0.015997
               8:          0.012405          0.011328          0.011677          0.012283
              15:          0.772434          0.075621          0.056711          0.079508
              16:          0.964635          0.105846          0.072811          0.104611
              31:         62.824051          2.295193          1.189505          2.351136
              32:         79.507060          1.169182          0.821779          1.275770
              63:           skipped         11.250225         10.982726         10.585744
              64:           skipped         10.013558         11.507645         12.665557



On a 2009 Mac Pro:
1-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000360          0.002269          0.004986          0.002476
               8:          0.000361          0.002468          0.005242          0.002696
              15:          0.000364          0.002255          0.005244          0.002471
              16:          0.000365          0.002506          0.005286          0.002727
              31:          0.000385          0.002380          0.005422          0.002588
              32:          0.000385          0.002531          0.005543          0.002737
              63:          0.000474          0.002407          0.005392          0.002637
              64:          0.000484          0.002602          0.005631          0.002823
             127:          0.000752          0.004122          0.007827          0.003966
             128:          0.000757          0.002763          0.005844          0.002958
             255:          0.004316          0.003258          0.006566          0.003324
             256:          0.004354          0.003180          0.006120          0.003245
             511:          0.011517          0.007158          0.009898          0.006238
             512:          0.011482          0.003873          0.006777          0.003820
            1023:          0.034105          0.009211          0.009468          0.008260
            1024:          0.034609          0.005504          0.008399          0.005080
            2047:          0.113620          0.028097          0.020662          0.021603
            2048:          0.112828          0.008403          0.010939          0.007331
            4095:          0.403373          0.023211          0.018767          0.020065
            4096:          0.403316          0.017550          0.017853          0.013651
            8191:          1.519329          8.454573          0.211436          7.212381
            8192:          1.519082          0.033148          0.030370          0.025905
           16383:          5.887481          0.317428          0.153344          0.237119
           16384:          5.888222          0.069379          0.065264          0.052847

2-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.000474          0.003470          0.006131          0.003503
               8:          0.000503          0.003565          0.006400          0.003586
              15:          0.002011          0.004481          0.007825          0.004496
              16:          0.002236          0.004744          0.007078          0.004680
              31:          0.027291          0.019433          0.014841          0.018034
              32:          0.029283          0.009244          0.010161          0.008964
              63:          0.445680          0.038171          0.026753          0.037404
              64:          0.460616          0.028128          0.029487          0.029149
             127:          7.003774          0.925921          0.282591          0.762671
             128:          7.063657          0.110838          0.104402          0.133523
             255:           skipped          0.804682          0.708849          0.869368
             256:           skipped          0.797800          0.721042          0.880848
             511:           skipped          3.643626          3.687562          4.584770
             512:           skipped          3.715215          4.893539          5.538462

3-dimensional arrays ('n' is the size of the image AND the kernel)
                n          convolve    convolve_fftnp     convolve_fftw    convolve_fftsp
               7:          0.004520          0.011519          0.009464          0.012335
               8:          0.006422          0.010294          0.010220          0.011711
              15:          0.329566          0.060978          0.045495          0.073692
              16:          0.405275          0.069999          0.040659          0.086114
              31:         24.935228          1.654920          0.710509          1.773879
              32:         27.524226          0.724053          0.543507          1.027568
              63:           skipped          8.982771         12.407683         16.900078
              64:           skipped          8.956070         11.934627         17.296447

"""
