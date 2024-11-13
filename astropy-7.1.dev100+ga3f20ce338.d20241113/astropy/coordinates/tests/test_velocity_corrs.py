import numpy as np
import pytest

from astropy import units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import Distance, EarthLocation, SkyCoord
from astropy.coordinates.sites import get_builtin_sites
from astropy.table import Table
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename


@pytest.fixture(scope="module")
def input_radecs():
    ras = []
    decs = []

    for dec in np.linspace(-85, 85, 15):
        nra = int(np.round(10 * np.cos(dec * u.deg)).value)
        ras1 = np.linspace(-180, 180 - 1e-6, nra)
        ras.extend(ras1)
        decs.extend([dec] * len(ras1))

    return SkyCoord(ra=ras, dec=decs, unit=u.deg)


@pytest.mark.parametrize("kind", ["heliocentric", "barycentric"])
def test_basic(kind):
    t0 = Time("2015-1-1")
    loc = get_builtin_sites()["example_site"]

    sc = SkyCoord(0, 0, unit=u.deg, obstime=t0, location=loc)
    rvc0 = sc.radial_velocity_correction(kind)

    assert rvc0.shape == ()
    assert rvc0.unit.is_equivalent(u.km / u.s)

    scs = SkyCoord(0, 0, unit=u.deg, obstime=t0 + np.arange(10) * u.day, location=loc)
    rvcs = scs.radial_velocity_correction(kind)
    assert rvcs.shape == (10,)
    assert rvcs.unit.is_equivalent(u.km / u.s)


test_input_time = Time(2457244.5, format="jd")
# test_input_loc = EarthLocation.of_site('Cerro Paranal')
# to avoid the network hit we just copy here what that yields
test_input_loc = EarthLocation.from_geodetic(
    lon=-70.403 * u.deg, lat=-24.6252 * u.deg, height=2635 * u.m
)


def test_helio_iraf(input_radecs):
    """
    Compare the heliocentric correction to the IRAF rvcorrect.
    `generate_IRAF_input` function is provided to show how the comparison data
    was produced

    def generate_IRAF_input(writefn=None):
        dt = test_input_time.utc.datetime

        coos = input_radecs  # `input_radecs` is implemented as pytest fixture

        lines = []
        for ra, dec in zip(coos.ra, coos.dec):
            rastr = Angle(ra).to_string(u.hour, sep=":")
            decstr = Angle(dec).to_string(u.deg, sep=":")

            lines.append(
                f"{dt.year} {dt.month} {dt.day} {dt.hour}:{dt.minute} {rastr} {decstr}"
            )
        if writefn:
            with open(writefn, "w") as f:
                for l in lines:
                    f.write(l)
        else:
            for l in lines:
                print(l)
        print("Run IRAF as:\nastutil\nrvcorrect f=<filename> observatory=Paranal")
    """
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
    for line in rvcorr_result.strip().split("\n")[5:]:
        vhs_iraf.append(float(line.split()[2]))
    vhs_iraf = vhs_iraf * u.km / u.s

    targets = SkyCoord(input_radecs, obstime=test_input_time, location=test_input_loc)
    vhs_astropy = targets.radial_velocity_correction("heliocentric")
    assert_quantity_allclose(vhs_astropy, vhs_iraf, atol=150 * u.m / u.s)


def test_barycorr(input_radecs):
    barycorr_bvcs = (
        np.loadtxt(get_pkg_data_filename("data/barycorr_bvcs.dat")) * u.m / u.s
    )

    # this tries the *other* way of calling radial_velocity_correction relative
    # to the IRAF tests
    bvcs_astropy = input_radecs.radial_velocity_correction(
        obstime=test_input_time, location=test_input_loc, kind="barycentric"
    )

    assert_quantity_allclose(bvcs_astropy, barycorr_bvcs, atol=10 * u.mm / u.s)


def test_rvcorr_multiple_obstimes_onskycoord():
    loc = EarthLocation(-2309223 * u.m, -3695529 * u.m, -4641767 * u.m)
    arrtime = Time("2005-03-21 00:00:00") + np.linspace(-1, 1, 10) * u.day

    sc = SkyCoord(1 * u.deg, 2 * u.deg, 100 * u.kpc, obstime=arrtime, location=loc)
    rvcbary_sc2 = sc.radial_velocity_correction(kind="barycentric")
    assert len(rvcbary_sc2) == 10

    # check the multiple-obstime and multi- mode
    sc = SkyCoord(
        ([1] * 10) * u.deg, 2 * u.deg, 100 * u.kpc, obstime=arrtime, location=loc
    )
    rvcbary_sc3 = sc.radial_velocity_correction(kind="barycentric")
    assert len(rvcbary_sc3) == 10


def test_invalid_argument_combos():
    loc = EarthLocation(-2309223 * u.m, -3695529 * u.m, -4641767 * u.m)
    time = Time("2005-03-21 00:00:00")
    timel = Time("2005-03-21 00:00:00", location=loc)

    scwattrs = SkyCoord(1 * u.deg, 2 * u.deg, obstime=time, location=loc)
    scwoattrs = SkyCoord(1 * u.deg, 2 * u.deg)

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
    sc = SkyCoord(
        10 * u.deg,
        20 * u.deg,
        distance=5 * u.pc,
        obstime=test_input_time,
        pm_ra_cosdec=0 * u.mas / u.yr,
        pm_dec=0 * u.mas / u.yr,
        radial_velocity=0 * u.km / u.s,
    )
    sc_novel = SkyCoord(
        10 * u.deg, 20 * u.deg, distance=5 * u.pc, obstime=test_input_time
    )
    corr = sc.radial_velocity_correction(
        obstime=test_input_time, location=test_input_loc
    )
    corr_novel = sc_novel.radial_velocity_correction(
        obstime=test_input_time, location=test_input_loc
    )
    assert_quantity_allclose(corr, corr_novel)


def test_barycorr_withvels(input_radecs):
    barycorr_bvcs = (
        np.loadtxt(get_pkg_data_filename("data/barycorr_bvcs_withvels.dat")) * u.m / u.s
    )
    bvcs_astropy = SkyCoord(
        input_radecs.ra,
        input_radecs.dec,
        pm_ra_cosdec=np.linspace(-1000, 1000, input_radecs.size) * u.mas / u.yr,
        pm_dec=np.linspace(0, 1000, input_radecs.size) * u.mas / u.yr,
        radial_velocity=np.linspace(0, 100, input_radecs.size) * u.km / u.s,
        distance=np.linspace(10, 100, input_radecs.size) * u.pc,
        obstime=test_input_time,
    ).radial_velocity_correction(obstime=test_input_time, location=test_input_loc)
    assert_quantity_allclose(bvcs_astropy, barycorr_bvcs, atol=10 * u.mm / u.s)


def test_warning_no_obstime_on_skycoord():
    c = SkyCoord(
        l=10 * u.degree,
        b=45 * u.degree,
        pm_l_cosb=34 * u.mas / u.yr,
        pm_b=-117 * u.mas / u.yr,
        distance=50 * u.pc,
        frame="galactic",
    )
    with pytest.warns(Warning):
        c.radial_velocity_correction("barycentric", test_input_time, test_input_loc)


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
        get_pkg_data_filename("coordinates/wright_eastmann_2014_tau_ceti.fits")
    )
    reduced_jds = wright_table["JD-2400000"]
    tempo2 = wright_table["TEMPO2"]
    barycorr = wright_table["BARYCORR"]

    # tau Ceti Hipparchos data
    tauCet = SkyCoord(
        "01 44 05.1275 -15 56 22.4006",
        unit=(u.hour, u.deg),
        pm_ra_cosdec=-1721.05 * u.mas / u.yr,
        pm_dec=854.16 * u.mas / u.yr,
        distance=Distance(parallax=273.96 * u.mas),
        radial_velocity=-16.597 * u.km / u.s,
        obstime=Time(48348.5625, format="mjd"),
    )
    # CTIO location as used in Wright & Eastmann
    xyz = u.Quantity([1814985.3, -5213916.8, -3187738.1], u.m)
    obs = EarthLocation(*xyz)
    times = Time(2400000, reduced_jds, format="jd")
    tempo2 = tempo2 * speed_of_light
    barycorr = barycorr * speed_of_light
    astropy = tauCet.radial_velocity_correction(location=obs, obstime=times)

    assert_quantity_allclose(astropy, tempo2, atol=5 * u.mm / u.s)
    assert_quantity_allclose(astropy, barycorr, atol=5 * u.mm / u.s)
