
import pytest

import numpy as np

from astropy.tests.helper import assert_quantity_allclose
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord, Angle, Distance
from astropy.coordinates.sites import get_builtin_sites
from astropy.utils.data import get_pkg_data_filename
from astropy.constants import c as speed_of_light
from astropy.table import Table


@pytest.mark.parametrize('kind', ['heliocentric', 'barycentric'])
def test_basic(kind):
    t0 = Time('2015-1-1')
    loc = get_builtin_sites()['example_site']

    sc = SkyCoord(0, 0, unit=u.deg, obstime=t0, location=loc)
    rvc0 = sc.radial_velocity_correction(kind)

    assert rvc0.shape == ()
    assert rvc0.unit.is_equivalent(u.km/u.s)

    scs = SkyCoord(0, 0, unit=u.deg, obstime=t0 + np.arange(10)*u.day,
                   location=loc)
    rvcs = scs.radial_velocity_correction(kind)
    assert rvcs.shape == (10,)
    assert rvcs.unit.is_equivalent(u.km/u.s)


test_input_time = Time(2457244.5, format='jd')
# test_input_loc = EarthLocation.of_site('Cerro Paranal')
# to avoid the network hit we just copy here what that yields
test_input_loc = EarthLocation.from_geodetic(lon=-70.403*u.deg,
                                             lat=-24.6252*u.deg,
                                             height=2635*u.m)


def test_helio_iraf():
    """
    Compare the heliocentric correction to the IRAF rvcorrect.
    `generate_IRAF_input` function is provided to show how the comparison data
    was produced

    """
    # this is based on running IRAF with the output of `generate_IRAF_input` below
    rvcorr_result = """
    # RVCORRECT: Observatory parameters for European Southern Observatory: Paranal
    #       latitude = -24:37.5
    #       longitude = 70:24.2
    #       altitude = 2635
    ##   HJD          VOBS   VHELIO     VLSR   VDIURNAL   VLUNAR  VANNUAL   VSOLAR
    2457244.50120     0.00   -10.36   -20.35     -0.034   -0.001  -10.325   -9.993
    2457244.50025     0.00   -14.20   -23.86     -0.115   -0.004  -14.085   -9.656
    2457244.50278     0.00    -2.29   -11.75      0.115    0.004   -2.413   -9.459
    2457244.50025     0.00   -14.20   -23.86     -0.115   -0.004  -14.085   -9.656
    2457244.49929     0.00   -17.41   -26.30     -0.192   -0.006  -17.214   -8.888
    2457244.50317     0.00   -17.19   -17.44      0.078    0.001  -17.269   -0.253
    2457244.50348     0.00     2.35    -6.21      0.192    0.006    2.156   -8.560
    2457244.49959     0.00     2.13   -15.06     -0.078   -0.000    2.211  -17.194
    2457244.49929     0.00   -17.41   -26.30     -0.192   -0.006  -17.214   -8.888
    2457244.49835     0.00   -19.84   -27.56     -0.259   -0.008  -19.573   -7.721
    2457244.50186     0.00   -24.47   -22.16     -0.038   -0.004  -24.433    2.313
    2457244.50470     0.00   -11.11    -8.57      0.221    0.005  -11.332    2.534
    2457244.50402     0.00     6.90    -0.38      0.259    0.008    6.629   -7.277
    2457244.50051     0.00    11.53    -5.78      0.038    0.004   11.489  -17.311
    2457244.49768     0.00    -1.84   -19.37     -0.221   -0.004   -1.612  -17.533
    2457244.49835     0.00   -19.84   -27.56     -0.259   -0.008  -19.573   -7.721
    2457244.49749     0.00   -21.38   -27.59     -0.315   -0.010  -21.056   -6.209
    2457244.50109     0.00   -27.69   -22.90     -0.096   -0.006  -27.584    4.785
    2457244.50457     0.00   -17.00    -9.30      0.196    0.003  -17.201    7.704
    2457244.50532     0.00     2.62     2.97      0.340    0.009    2.276    0.349
    2457244.50277     0.00    16.42     4.67      0.228    0.009   16.178  -11.741
    2457244.49884     0.00    13.98    -5.48     -0.056    0.002   14.039  -19.463
    2457244.49649     0.00    -2.84   -19.84     -0.297   -0.007   -2.533  -17.000
    2457244.49749     0.00   -21.38   -27.59     -0.315   -0.010  -21.056   -6.209
    2457244.49675     0.00   -21.97   -26.39     -0.357   -0.011  -21.598   -4.419
    2457244.50025     0.00   -29.30   -22.47     -0.149   -0.008  -29.146    6.831
    2457244.50398     0.00   -21.55    -9.88      0.146    0.001  -21.700   11.670
    2457244.50577     0.00    -3.26     4.00      0.356    0.009   -3.623    7.263
    2457244.50456     0.00    14.87    11.06      0.357    0.011   14.497   -3.808
    2457244.50106     0.00    22.20     7.14      0.149    0.008   22.045  -15.058
    2457244.49732     0.00    14.45    -5.44     -0.146   -0.001   14.600  -19.897
    2457244.49554     0.00    -3.84   -19.33     -0.356   -0.008   -3.478  -15.491
    2457244.49675     0.00   -21.97   -26.39     -0.357   -0.011  -21.598   -4.419
    2457244.49615     0.00   -21.57   -24.00     -0.383   -0.012  -21.172   -2.432
    2457244.49942     0.00   -29.36   -20.83     -0.193   -0.009  -29.157    8.527
    2457244.50312     0.00   -24.26    -9.75      0.088   -0.001  -24.348   14.511
    2457244.50552     0.00    -8.66     4.06      0.327    0.007   -8.996   12.721
    2457244.50549     0.00    10.14    14.13      0.413    0.012    9.715    3.994
    2457244.50305     0.00    23.35    15.76      0.306    0.011   23.031   -7.586
    2457244.49933     0.00    24.78     8.18      0.056    0.006   24.721  -16.601
    2457244.49609     0.00    13.77    -5.06     -0.221   -0.003   13.994  -18.832
    2457244.49483     0.00    -4.53   -17.77     -0.394   -0.010   -4.131  -13.237
    2457244.49615     0.00   -21.57   -24.00     -0.383   -0.012  -21.172   -2.432
    2457244.49572     0.00   -20.20   -20.54     -0.392   -0.013  -19.799   -0.335
    2457244.49907     0.00   -28.17   -17.30     -0.197   -0.009  -27.966   10.874
    2457244.50285     0.00   -22.96    -5.96      0.090   -0.001  -23.048   16.995
    2457244.50531     0.00    -7.00     8.16      0.335    0.007   -7.345   15.164
    2457244.50528     0.00    12.23    18.47      0.423    0.012   11.795    6.238
    2457244.50278     0.00    25.74    20.13      0.313    0.012   25.416   -5.607
    2457244.49898     0.00    27.21    12.38      0.057    0.006   27.144  -14.829
    2457244.49566     0.00    15.94    -1.17     -0.226   -0.003   16.172  -17.111
    2457244.49437     0.00    -2.78   -14.17     -0.403   -0.010   -2.368  -11.387
    2457244.49572     0.00   -20.20   -20.54     -0.392   -0.013  -19.799   -0.335
    2457244.49548     0.00   -17.94   -16.16     -0.383   -0.012  -17.541    1.776
    2457244.49875     0.00   -25.73   -12.99     -0.193   -0.009  -25.525   12.734
    2457244.50246     0.00   -20.63    -1.91      0.088   -0.001  -20.716   18.719
    2457244.50485     0.00    -5.03    11.90      0.327    0.007   -5.365   16.928
    2457244.50482     0.00    13.77    21.97      0.413    0.012   13.347    8.202
    2457244.50238     0.00    26.98    23.60      0.306    0.011   26.663   -3.378
    2457244.49867     0.00    28.41    16.02      0.056    0.005   28.353  -12.393
    2457244.49542     0.00    17.40     2.78     -0.221   -0.003   17.625  -14.625
    2457244.49416     0.00    -0.90    -9.93     -0.394   -0.010   -0.499   -9.029
    2457244.49548     0.00   -17.94   -16.16     -0.383   -0.012  -17.541    1.776
    2457244.49544     0.00   -14.87   -11.06     -0.357   -0.011  -14.497    3.808
    2457244.49894     0.00   -22.20    -7.14     -0.149   -0.008  -22.045   15.058
    2457244.50268     0.00   -14.45     5.44      0.146    0.001  -14.600   19.897
    2457244.50446     0.00     3.84    19.33      0.356    0.008    3.478   15.491
    2457244.50325     0.00    21.97    26.39      0.357    0.011   21.598    4.419
    2457244.49975     0.00    29.30    22.47      0.149    0.008   29.146   -6.831
    2457244.49602     0.00    21.55     9.88     -0.146   -0.001   21.700  -11.670
    2457244.49423     0.00     3.26    -4.00     -0.356   -0.009    3.623   -7.263
    2457244.49544     0.00   -14.87   -11.06     -0.357   -0.011  -14.497    3.808
    2457244.49561     0.00   -11.13    -5.46     -0.315   -0.010  -10.805    5.670
    2457244.49921     0.00   -17.43    -0.77     -0.096   -0.006  -17.333   16.664
    2457244.50269     0.00    -6.75    12.83      0.196    0.003   -6.949   19.583
    2457244.50344     0.00    12.88    25.10      0.340    0.009   12.527   12.227
    2457244.50089     0.00    26.67    26.80      0.228    0.009   26.430    0.137
    2457244.49696     0.00    24.24    16.65     -0.056    0.002   24.290   -7.584
    2457244.49461     0.00     7.42     2.29     -0.297   -0.007    7.719   -5.122
    2457244.49561     0.00   -11.13    -5.46     -0.315   -0.010  -10.805    5.670
    2457244.49598     0.00    -6.90     0.38     -0.259   -0.008   -6.629    7.277
    2457244.49949     0.00   -11.53     5.78     -0.038   -0.004  -11.489   17.311
    2457244.50232     0.00     1.84    19.37      0.221    0.004    1.612   17.533
    2457244.50165     0.00    19.84    27.56      0.259    0.008   19.573    7.721
    2457244.49814     0.00    24.47    22.16      0.038    0.004   24.433   -2.313
    2457244.49530     0.00    11.11     8.57     -0.221   -0.005   11.332   -2.534
    2457244.49598     0.00    -6.90     0.38     -0.259   -0.008   -6.629    7.277
    2457244.49652     0.00    -2.35     6.21     -0.192   -0.006   -2.156    8.560
    2457244.50041     0.00    -2.13    15.06      0.078    0.000   -2.211   17.194
    2457244.50071     0.00    17.41    26.30      0.192    0.006   17.214    8.888
    2457244.49683     0.00    17.19    17.44     -0.078   -0.001   17.269    0.253
    2457244.49652     0.00    -2.35     6.21     -0.192   -0.006   -2.156    8.560
    2457244.49722     0.00     2.29    11.75     -0.115   -0.004    2.413    9.459
    2457244.49975     0.00    14.20    23.86      0.115    0.004   14.085    9.656
    2457244.49722     0.00     2.29    11.75     -0.115   -0.004    2.413    9.459
    2457244.49805     0.00     6.84    16.77     -0.034   -0.001    6.874    9.935
    """
    vhs_iraf = []
    for line in rvcorr_result.strip().split('\n'):
        if not line.strip().startswith('#'):
            vhs_iraf.append(float(line.split()[2]))
    vhs_iraf = vhs_iraf*u.km/u.s

    targets = SkyCoord(_get_test_input_radecs(), obstime=test_input_time,
                       location=test_input_loc)
    vhs_astropy = targets.radial_velocity_correction('heliocentric')
    assert_quantity_allclose(vhs_astropy, vhs_iraf, atol=150*u.m/u.s)
    return vhs_astropy, vhs_iraf  # for interactively examination


def generate_IRAF_input(writefn=None):
    dt = test_input_time.utc.datetime

    coos = _get_test_input_radecs()

    lines = []
    for ra, dec in zip(coos.ra, coos.dec):
        rastr = Angle(ra).to_string(u.hour, sep=':')
        decstr = Angle(dec).to_string(u.deg, sep=':')

        msg = '{yr} {mo} {day} {uth}:{utmin} {ra} {dec}'
        lines.append(msg.format(yr=dt.year, mo=dt.month, day=dt.day,
                                uth=dt.hour, utmin=dt.minute,
                                ra=rastr, dec=decstr))
    if writefn:
        with open(writefn, 'w') as f:
            for l in lines:
                f.write(l)
    else:
        for l in lines:
            print(l)
    print('Run IRAF as:\nastutil\nrvcorrect f=<filename> observatory=Paranal')


def _get_test_input_radecs():
    ras = []
    decs = []

    for dec in np.linspace(-85, 85, 15):
        nra = int(np.round(10*np.cos(dec*u.deg)).value)
        ras1 = np.linspace(-180, 180-1e-6, nra)
        ras.extend(ras1)
        decs.extend([dec]*len(ras1))

    return SkyCoord(ra=ras, dec=decs, unit=u.deg)


def test_barycorr():
    # this is the result of calling _get_barycorr_bvcs
    barycorr_bvcs = u.Quantity([
       -10335.93326096, -14198.47605491, -2237.60012494, -14198.47595363,
       -17425.46512587, -17131.70901174, 2424.37095076, 2130.61519166,
       -17425.46495779, -19872.50026998, -24442.37091097, -11017.08975893,
         6978.0622355, 11547.93333743, -1877.34772637, -19872.50004258,
       -21430.08240017, -27669.14280689, -16917.08506807, 2729.57222968,
        16476.49569232, 13971.97171764, -2898.04250914, -21430.08212368,
       -22028.51337105, -29301.92349394, -21481.13036199, -3147.44828909,
        14959.50065514, 22232.91155425, 14412.11903105, -3921.56359768,
       -22028.51305781, -21641.01479409, -29373.0512649, -24205.90521765,
        -8557.34138828, 10250.50350732, 23417.2299926, 24781.98057941,
        13706.17339044, -4627.70005932, -21641.01445812, -20284.92627505,
       -28193.91696959, -22908.51624166, -6901.82132125, 12336.45758056,
        25804.51614607, 27200.50029664, 15871.21385688, -2882.24738355,
       -20284.9259314, -18020.92947805, -25752.96564978, -20585.81957567,
        -4937.25573801, 13870.58916957, 27037.31568441, 28402.06636994,
        17326.25977035, -1007.62209045, -18020.92914212, -14950.33284575,
       -22223.74260839, -14402.94943965, 3930.73265119, 22037.68163353,
        29311.09265126, 21490.30070307, 3156.62229843, -14950.33253252,
       -11210.53846867, -17449.59867676, -6697.54090389, 12949.11642965,
        26696.03999586, 24191.5164355, 7321.50355488, -11210.53819218,
        -6968.89359681, -11538.76423011, 1886.51695238, 19881.66902396,
        24451.54039956, 11026.26000765, -6968.89336945, -2415.20201758,
        -2121.44599781, 17434.63406085, 17140.87871753, -2415.2018495,
         2246.76923076, 14207.64513054, 2246.76933194, 6808.40787728],
         u.m/u.s)

    # this tries the *other* way of calling radial_velocity_correction relative
    # to the IRAF tests
    targets = _get_test_input_radecs()
    bvcs_astropy = targets.radial_velocity_correction(obstime=test_input_time,
                                                      location=test_input_loc,
                                                      kind='barycentric')

    assert_quantity_allclose(bvcs_astropy, barycorr_bvcs, atol=10*u.mm/u.s)
    return bvcs_astropy, barycorr_bvcs  # for interactively examination


def _get_barycorr_bvcs(coos, loc, injupyter=False):
    """
    Gets the barycentric correction of the test data from the
    http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.html web site.
    Requires the https://github.com/tronsgaard/barycorr python interface to that
    site.

    Provided to reproduce the test data above, but not required to actually run
    the tests.
    """
    import barycorr
    from astropy.utils.console import ProgressBar

    bvcs = []
    for ra, dec in ProgressBar(list(zip(coos.ra.deg, coos.dec.deg)),
                               ipython_widget=injupyter):
        res = barycorr.bvc(test_input_time.utc.jd, ra, dec,
                           lat=loc.geodetic[1].deg,
                           lon=loc.geodetic[0].deg,
                           elevation=loc.geodetic[2].to(u.m).value)
        bvcs.append(res)
    return bvcs*u.m/u.s


def test_rvcorr_multiple_obstimes_onskycoord():
    loc = EarthLocation(-2309223 * u.m, -3695529 * u.m, -4641767 * u.m)
    arrtime = Time('2005-03-21 00:00:00') + np.linspace(-1, 1, 10)*u.day

    sc = SkyCoord(1*u.deg, 2*u.deg, 100*u.kpc, obstime=arrtime, location=loc)
    rvcbary_sc2 = sc.radial_velocity_correction(kind='barycentric')
    assert len(rvcbary_sc2) == 10

    # check the multiple-obstime and multi- mode
    sc = SkyCoord(([1]*10)*u.deg, 2*u.deg, 100*u.kpc,
                  obstime=arrtime, location=loc)
    rvcbary_sc3 = sc.radial_velocity_correction(kind='barycentric')
    assert len(rvcbary_sc3) == 10


def test_invalid_argument_combos():
    loc = EarthLocation(-2309223 * u.m, -3695529 * u.m, -4641767 * u.m)
    time = Time('2005-03-21 00:00:00')
    timel = Time('2005-03-21 00:00:00', location=loc)

    scwattrs = SkyCoord(1*u.deg, 2*u.deg, obstime=time, location=loc)
    scwoattrs = SkyCoord(1*u.deg, 2*u.deg)

    scwattrs.radial_velocity_correction()
    with pytest.raises(ValueError):
        scwattrs.radial_velocity_correction(obstime=time, location=loc)
    with pytest.raises(TypeError):
        scwoattrs.radial_velocity_correction(obstime=time)

    scwoattrs.radial_velocity_correction(obstime=time, location=loc)
    with pytest.raises(TypeError):
        scwoattrs.radial_velocity_correction()

    with pytest.raises(ValueError):
        scwattrs.radial_velocity_correction(timel)


def test_regression_9645():
    sc = SkyCoord(10*u.deg, 20*u.deg, distance=5*u.pc, obstime=test_input_time,
                  pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr, radial_velocity=0*u.km/u.s)
    sc_novel = SkyCoord(10*u.deg, 20*u.deg, distance=5*u.pc, obstime=test_input_time)
    corr = sc.radial_velocity_correction(obstime=test_input_time, location=test_input_loc)
    corr_novel = sc_novel.radial_velocity_correction(obstime=test_input_time, location=test_input_loc)
    assert_quantity_allclose(corr, corr_novel)


def test_barycorr_withvels():
    # this is the result of calling _get_barycorr_bvcs_withvels
    barycorr_bvcs = u.Quantity(
        [-10335.94926581, -14198.49117304,  -2237.58656335,
         -14198.49078575, -17425.47883864, -17131.72711182,
         2424.38466675,   2130.62819093, -17425.47834604,
         -19872.51254565, -24442.39064348, -11017.0964353,
         6978.07515501,  11547.94831175,  -1877.34560543,
         -19872.51188308, -21430.0931411, -27669.15919972,
         -16917.09482078,   2729.57757823,  16476.5087925,
         13971.97955641,  -2898.04451551, -21430.09220144,
         -22028.52224227, -29301.93613248, -21481.14015151,
         -3147.44852058,  14959.50849997,  22232.91906264,
         14412.12044201,  -3921.56783473, -22028.52088749,
         -21641.02117064, -29373.05982792, -24205.91319258,
         -8557.34473049,  10250.50560918,  23417.23357219,
         24781.98113432,  13706.17025059,  -4627.70468688,
         -21641.01928189, -20284.92926795, -28193.92117514,
         -22908.52127321,  -6901.82512637,  12336.45557256,
         25804.5137786,  27200.49576347,  15871.20847332,
         -2882.25080211, -20284.92696256, -18020.92824383,
         -25752.96528309, -20585.82211189,  -4937.26088706,
         13870.58217495,  27037.30698639,  28402.0571686,
         17326.25314311,  -1007.62313006, -18020.92552769,
         -14950.32653444, -22223.73793506, -14402.95155047,
         3930.72325162,  22037.66749783,  29311.07826101,
         21490.29193529,   3156.62360741, -14950.32373745,
         -11210.52665171, -17449.59068509,  -6697.54579192,
         12949.09948082,  26696.01956077,  24191.50403015,
         7321.50684816, -11210.52389393,  -6968.87610888,
         -11538.7547047,   1886.50525065,  19881.64366561,
         24451.52197666,  11026.26396455,  -6968.87351156,
         -2415.17899385,  -2121.44598968,  17434.60465075,
         17140.87204017,  -2415.1771038,   2246.79688215,
         14207.61339552,   2246.79790276,   6808.43888253], u.m/u.s)
    coos = _get_test_input_radecvels()
    bvcs_astropy = coos.radial_velocity_correction(obstime=test_input_time,
                                                   location=test_input_loc)
    assert_quantity_allclose(bvcs_astropy, barycorr_bvcs, atol=10*u.mm/u.s)
    return bvcs_astropy, barycorr_bvcs  # for interactively examination


def _get_test_input_radecvels():
    coos = _get_test_input_radecs()
    ras = coos.ra
    decs = coos.dec
    pmra = np.linspace(-1000, 1000, coos.size)*u.mas/u.yr
    pmdec = np.linspace(0, 1000, coos.size)*u.mas/u.yr
    rvs = np.linspace(0, 100, coos.size)*u.km/u.s
    distance = np.linspace(10, 100, coos.size)*u.pc
    return SkyCoord(ras, decs, pm_ra_cosdec=pmra, pm_dec=pmdec,
                    radial_velocity=rvs, distance=distance,
                    obstime=test_input_time)


def _get_barycorr_bvcs_withvels(coos, loc, injupyter=False):
    """
    Gets the barycentric correction of the test data from the
    http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.html web site.
    Requires the https://github.com/tronsgaard/barycorr python interface to that
    site.

    Provided to reproduce the test data above, but not required to actually run
    the tests.
    """
    import barycorr
    from astropy.utils.console import ProgressBar

    bvcs = []
    for coo in ProgressBar(coos, ipython_widget=injupyter):
        res = barycorr.bvc(test_input_time.utc.jd,
                           coo.ra.deg, coo.dec.deg,
                           lat=loc.geodetic[1].deg,
                           lon=loc.geodetic[0].deg,
                           pmra=coo.pm_ra_cosdec.to_value(u.mas/u.yr),
                           pmdec=coo.pm_dec.to_value(u.mas/u.yr),
                           parallax=coo.distance.to_value(u.mas, equivalencies=u.parallax()),
                           rv=coo.radial_velocity.to_value(u.m/u.s),
                           epoch=test_input_time.utc.jd,
                           elevation=loc.geodetic[2].to(u.m).value)
        bvcs.append(res)
    return bvcs*u.m/u.s


def test_warning_no_obstime_on_skycoord():
    c = SkyCoord(l=10*u.degree, b=45*u.degree,
                 pm_l_cosb=34*u.mas/u.yr, pm_b=-117*u.mas/u.yr,
                 distance=50*u.pc, frame='galactic')
    with pytest.warns(Warning):
        c.radial_velocity_correction('barycentric', test_input_time,
                                     test_input_loc)


@pytest.mark.remote_data
def test_regression_10094():
    """
    Make sure that when we include the proper motion and radial velocity of
    a SkyCoord, our velocity corrections remain close to TEMPO2.

    We check that tau Ceti is within 5mm/s
    """
    # Wright & Eastman (2014) Table2
    # Corrections for tau Ceti

    wright_table = Table.read(
        get_pkg_data_filename('coordinates/wright_eastmann_2014_tau_ceti.fits')
    )
    reduced_jds = wright_table['JD-2400000']
    tempo2 = wright_table['TEMPO2']
    barycorr = wright_table['BARYCORR']

    # tau Ceti Hipparchos data
    tauCet = SkyCoord('01 44 05.1275 -15 56 22.4006',
                      unit=(u.hour, u.deg),
                      pm_ra_cosdec=-1721.05*u.mas/u.yr,
                      pm_dec=854.16*u.mas/u.yr,
                      distance=Distance(parallax=273.96*u.mas),
                      radial_velocity=-16.597*u.km/u.s,
                      obstime=Time(48348.5625, format='mjd'))
    # CTIO location as used in Wright & Eastmann
    xyz = u.Quantity([1814985.3, -5213916.8, -3187738.1], u.m)
    obs = EarthLocation(*xyz)
    times = Time(2400000, reduced_jds, format='jd')
    tempo2 = tempo2 * speed_of_light
    barycorr = barycorr * speed_of_light
    astropy = tauCet.radial_velocity_correction(location=obs, obstime=times)

    assert_quantity_allclose(astropy, tempo2, atol=5*u.mm/u.s)
    assert_quantity_allclose(astropy, barycorr, atol=5*u.mm/u.s)
