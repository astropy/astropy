# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.tests.helper import assert_quantity_allclose
from astropy.timeseries.periodograms.bls import BoxLeastSquares
from astropy.timeseries.periodograms.lombscargle.core import has_units


def assert_allclose_blsresults(blsresult, other, **kwargs):
    """Assert that another BoxLeastSquaresResults object is consistent

    This method loops over all attributes and compares the values using
    :func:`~astropy.tests.helper.assert_quantity_allclose` function.

    Parameters
    ----------
    other : BoxLeastSquaresResults
        The other results object to compare.

    """
    for k, v in blsresult.items():
        if k not in other:
            raise AssertionError(f"missing key '{k}'")
        if k == "objective":
            assert v == other[k], (
                f"Mismatched objectives. Expected '{v}', got '{other[k]}'"
            )
            continue
        assert_quantity_allclose(v, other[k], **kwargs)


# NOTE: PR 10644 replaced deprecated usage of RandomState but could not
# find a new seed that did not cause test failure, resorted to hardcoding.
@pytest.fixture
def data():
    t = np.array([
        6.96469186, 2.86139335, 2.26851454, 5.51314769, 7.1946897,
        4.2310646, 9.80764198, 6.84829739, 4.80931901, 3.92117518,
        3.43178016, 7.29049707, 4.38572245, 0.59677897, 3.98044255,
        7.37995406, 1.8249173, 1.75451756, 5.31551374, 5.31827587,
        6.34400959, 8.49431794, 7.24455325, 6.11023511, 7.22443383,
        3.22958914, 3.61788656, 2.28263231, 2.93714046, 6.30976124,
        0.9210494, 4.33701173, 4.30862763, 4.93685098, 4.2583029,
        3.12261223, 4.26351307, 8.93389163, 9.44160018, 5.01836676,
        6.23952952, 1.15618395, 3.17285482, 4.14826212, 8.66309158,
        2.50455365, 4.83034264, 9.85559786, 5.19485119, 6.12894526,
        1.20628666, 8.26340801, 6.03060128, 5.45068006, 3.42763834,
        3.04120789, 4.17022211, 6.81300766, 8.75456842, 5.10422337,
        6.69313783, 5.85936553, 6.24903502, 6.74689051, 8.42342438,
        0.83194988, 7.63682841, 2.43666375, 1.94222961, 5.72456957,
        0.95712517, 8.85326826, 6.27248972, 7.23416358, 0.16129207,
        5.94431879, 5.56785192, 1.58959644, 1.53070515, 6.95529529,
        3.18766426, 6.91970296, 5.5438325, 3.88950574, 9.2513249,
        8.41669997, 3.57397567, 0.43591464, 3.04768073, 3.98185682,
        7.0495883, 9.95358482, 3.55914866, 7.62547814, 5.93176917,
        6.91701799, 1.51127452, 3.98876293, 2.40855898, 3.43456014,
        5.13128154, 6.6662455, 1.05908485, 1.30894951, 3.21980606,
        6.61564337, 8.46506225, 5.53257345, 8.54452488, 3.84837811,
        3.16787897, 3.54264676, 1.71081829, 8.29112635, 3.38670846,
        5.52370075, 5.78551468, 5.21533059, 0.02688065, 9.88345419,
        9.05341576, 2.07635861, 2.92489413, 5.20010153, 9.01911373,
        9.83630885, 2.57542064, 5.64359043, 8.06968684, 3.94370054,
        7.31073036, 1.61069014, 6.00698568, 8.65864458, 9.83521609,
        0.7936579, 4.28347275, 2.0454286, 4.50636491, 5.47763573,
        0.9332671, 2.96860775, 9.2758424, 5.69003731, 4.57411998,
        7.53525991, 7.41862152, 0.48579033, 7.08697395, 8.39243348,
        1.65937884, 7.80997938, 2.86536617, 3.06469753, 6.65261465,
        1.11392172, 6.64872449, 8.87856793, 6.96311268, 4.40327877,
        4.38214384, 7.65096095, 5.65642001, 0.84904163, 5.82671088,
        8.14843703, 3.37066383, 9.2757658, 7.50717, 5.74063825,
        7.51643989, 0.79148961, 8.59389076, 8.21504113, 9.0987166,
        1.28631198, 0.81780087, 1.38415573, 3.9937871, 4.24306861,
        5.62218379, 1.2224355, 2.01399501, 8.11644348, 4.67987574,
        8.07938209, 0.07426379, 5.51592726, 9.31932148, 5.82175459,
        2.06095727, 7.17757562, 3.7898585, 6.68383947, 0.29319723,
        6.35900359, 0.32197935, 7.44780655, 4.72913002, 1.21754355,
        5.42635926, 0.66774443, 6.53364871, 9.96086327, 7.69397337,
        5.73774114, 1.02635259, 6.99834075, 6.61167867, 0.49097131,
        7.92299302, 5.18716591, 4.25867694, 7.88187174, 4.11569223,
        4.81026276, 1.81628843, 3.213189, 8.45532997, 1.86903749,
        4.17291061, 9.89034507, 2.36599812, 9.16832333, 9.18397468,
        0.91296342, 4.63652725, 5.02216335, 3.1366895, 0.47339537,
        2.41685637, 0.95529642, 2.38249906, 8.07791086, 8.94978288,
        0.43222892, 3.01946836, 9.80582199, 5.39504823, 6.26309362,
        0.05545408, 4.84909443, 9.88328535, 3.75185527, 0.97038159,
        4.61908762, 9.63004466, 3.41830614, 7.98922733, 7.98846331,
        2.08248297, 4.43367702, 7.15601275, 4.10519785, 1.91006955,
        9.67494307, 6.50750366, 8.65459852, 2.52423578e-01, 2.66905815,
        5.02071100, 6.74486351e-01, 9.93033261, 2.36462396, 3.74292182,
        2.14011915, 1.05445866, 2.32479786, 3.00610136, 6.34442268,
        2.81234781, 3.62276761, 5.94284372e-02, 3.65719126, 5.33885982,
        1.62015837, 5.97433108, 2.93152469, 6.32050495, 2.61966053e-01,
        8.87593460, 1.61186304e-01, 1.26958031, 7.77162462, 4.58952322e-01,
        7.10998694, 9.71046141, 8.71682933, 7.10161651, 9.58509743,
        4.29813338, 8.72878914, 3.55957668, 9.29763653, 1.48777656,
        9.40029015, 8.32716197, 8.46054838, 1.23923010, 5.96486898,
        1.63924809e-01, 7.21184366, 7.73751413e-02, 8.48222774e-01, 2.25498410,
        8.75124534, 3.63576318, 5.39959935, 5.68103214, 2.25463360,
        5.72146768, 6.60951795, 2.98245393, 4.18626859, 4.53088925,
        9.32350662, 5.87493747, 9.48252372, 5.56034754, 5.00561421,
        3.53221097e-02, 4.80889044, 9.27454999, 1.98365689, 5.20911344e-01,
        4.06778893, 3.72396481, 8.57153058, 2.66111156e-01, 9.20149230,
        6.80902999, 9.04225994, 6.07529071, 8.11953312, 3.35543874,
        3.49566228, 3.89874230, 7.54797082, 3.69291174, 2.42219806,
        9.37668357, 9.08011084, 3.48797316, 6.34638070, 2.73842212,
        2.06115129, 3.36339529, 3.27099893, 8.82276101, 8.22303815,
        7.09623229, 9.59345225, 4.22543353, 2.45033039, 1.17398437,
        3.01053358, 1.45263734, 9.21860974e-01, 6.02932197, 3.64187450,
        5.64570343, 1.91335721, 6.76905860, 2.15505447, 2.78023594,
        7.41760422, 5.59737896, 3.34836413, 5.42988783, 6.93984703,
        9.12132121, 5.80713213, 2.32686379, 7.46697631, 7.77769018,
        2.00401315, 8.20574220, 4.64934855, 7.79766662, 2.37478220,
        3.32580270, 9.53697119, 6.57815073, 7.72877831, 6.88374343,
        2.04304118, 4.70688748, 8.08963873, 6.75035127, 6.02788565e-02,
        8.74077427e-01, 3.46794720, 9.44365540, 4.91190481, 2.70176267,
        3.60423719, 2.10652628, 4.21200057, 2.18035440, 8.45752507,
        4.56270599, 2.79802018, 9.32891648, 3.14351354, 9.09714662,
        4.34180910e-01, 7.07115060, 4.83889039, 4.44221061, 3.63233444e-01,
        4.06831905e-01, 3.32753617, 9.47119540, 6.17659977, 3.68874842,
        6.11977039, 2.06131536, 1.65066443, 3.61817266, 8.63353352,
        5.09401727, 2.96901516, 9.50251625, 8.15966090, 3.22973943,
        9.72098245, 9.87351098, 4.08660134, 6.55923103, 4.05653198,
        2.57348106, 8.26526760e-01, 2.63610346, 2.71479854, 3.98639080,
        1.84886031, 9.53818403, 1.02879885,  6.25208533, 4.41697388,
        4.23518049, 3.71991783, 8.68314710, 2.80476981, 2.05761574e-01,
        9.18097016, 8.64480278, 2.76901790, 5.23487548, 1.09088197,
        9.34270688e-01, 8.37466108, 4.10265718, 6.61716540, 9.43200558,
        2.45130592, 1.31598313e-01, 2.41484058e-01, 7.09385692, 9.24551885,
        4.67330273, 3.75109148, 5.42860425, 8.58916838, 6.52153874,
        2.32979897, 7.74580205, 1.34613497, 1.65559971, 6.12682283,
        2.38783406, 7.04778548, 3.49518527, 2.77423960, 9.98918406,
        4.06161246e-01, 6.45822522, 3.86995850e-01, 7.60210258, 2.30089957,
        8.98318671e-01, 6.48449712, 7.32601217, 6.78095315, 5.19009471e-01,
        2.94306946, 4.51088346, 2.87103290, 8.10513456, 1.31115105,
        6.12179362, 9.88214944, 9.02556539, 2.22157062, 8.18876137e-04,
        9.80597342, 8.82712985, 9.19472466, 4.15503551, 7.44615462])
    y = np.ones_like(t)
    dy = np.array([
        0.00606416, 0.00696152, 0.00925774, 0.00563806, 0.00946933,
        0.00748254, 0.00713048, 0.00652823, 0.00958424, 0.00758812,
        0.00902013, 0.00928826, 0.00961191, 0.0065169, 0.00669905,
        0.00797537, 0.00720662, 0.00966421, 0.00698782, 0.00738889,
        0.00808593, 0.0070237, 0.00996239, 0.00549426, 0.00610302,
        0.00661328, 0.00573861, 0.0064211, 0.00889623, 0.00761446,
        0.00516977, 0.00991311, 0.00808003, 0.0052947, 0.00830584,
        0.00689185, 0.00567837, 0.00781832, 0.0086354, 0.00835563,
        0.00623757, 0.00762433, 0.00768832, 0.00858402, 0.00679934,
        0.00898866, 0.00813961, 0.00519166, 0.0077324, 0.00930956,
        0.00783787, 0.00587914, 0.00755188, 0.00878473, 0.00555053,
        0.0090855, 0.00583741, 0.00767038, 0.00692872, 0.00624312,
        0.00823716, 0.00518696, 0.00880023, 0.0076347, 0.00937886,
        0.00760359, 0.00517517, 0.005718, 0.00897802, 0.00745988,
        0.0072094, 0.00659217, 0.00642275, 0.00982943, 0.00716485,
        0.00942002, 0.00824082, 0.00929214, 0.00926225, 0.00978156,
        0.00848971, 0.00902698, 0.00866564, 0.00802613, 0.00858677,
        0.00857875, 0.00520454, 0.00758055, 0.00896326, 0.00621481,
        0.00732574, 0.00717493, 0.00701394, 0.0056092, 0.00762856,
        0.00723124, 0.00831696, 0.00774707, 0.00513771, 0.00515959,
        0.0085068, 0.00853791, 0.0097997, 0.00938352, 0.0073403,
        0.00812953, 0.00728591, 0.00611473, 0.00688338, 0.00551942,
        0.00833264, 0.00596015, 0.00737734, 0.00983718, 0.00515834,
        0.00575865, 0.0064929, 0.00970903, 0.00954421, 0.00581,
        0.00990559, 0.00875374, 0.00769989, 0.00965851, 0.00940304,
        0.00695658, 0.00828172, 0.00823693, 0.00663484, 0.00589695,
        0.00733405, 0.00631641, 0.00677533, 0.00977072, 0.00730569,
        0.00842446, 0.00668115, 0.00997931, 0.00829384, 0.00598005,
        0.00549092, 0.0097159, 0.00972389, 0.00810664, 0.00508496,
        0.00612767, 0.00900638, 0.0093773, 0.00726995, 0.0068276,
        0.00637113, 0.00558485, 0.00557872, 0.00976301, 0.00904313,
        0.0058239, 0.00603525, 0.00827776, 0.00882332, 0.00905157,
        0.00581669, 0.00992064, 0.00613901, 0.00794708, 0.00793808,
        0.00983681, 0.00828834, 0.00792452, 0.00759386, 0.00882329,
        0.00553028, 0.00501046, 0.00976244, 0.00749329, 0.00664168,
        0.00684027, 0.00901922, 0.00691185, 0.00885085, 0.00720231,
        0.00922039, 0.00538102, 0.00740564, 0.00733425, 0.00632164,
        0.00971807, 0.00952514, 0.00721798, 0.0054858, 0.00603392,
        0.00635746, 0.0074211, 0.00669189, 0.00887068, 0.00738013,
        0.00935185, 0.00997891, 0.00609918, 0.00805836, 0.00923751,
        0.00972618, 0.00645043, 0.00863521, 0.00507508, 0.00939571,
        0.00531969, 0.00866698, 0.00997305, 0.00750595, 0.00604667,
        0.00797322, 0.00812075, 0.00834036, 0.00586306, 0.00949356,
        0.00810496, 0.00521784, 0.00842021, 0.00598042, 0.0051367,
        0.00775477, 0.00906657, 0.00929971, 0.0055176, 0.00831521,
        0.00855038, 0.00647258, 0.00985682, 0.00639344, 0.00534991,
        0.0075964, 0.00847157, 0.0062233, 0.00669291, 0.00781814,
        0.00943339, 0.00873663, 0.00604796, 0.00625889, 0.0076194,
        0.00884479, 0.00809381, 0.00750662, 0.00798563, 0.0087803,
        0.0076854, 0.00948876, 0.00973534, 0.00957677, 0.00877259,
        0.00623161, 0.00692636, 0.0064, 0.0082883, 0.00662111,
        0.00877196, 0.00556755, 0.00887682, 0.00792951, 0.00917694,
        0.00715438, 0.00812482, 0.00777206, 0.00987836, 0.00877737,
        0.00772407, 0.00587016, 0.00952057, 0.00602919, 0.00825022,
        0.00968236, 0.0061179, 0.00612962, 0.00925909, 0.00913828,
        0.00675852, 0.00632548, 0.00563694, 0.00993968, 0.00917672,
        0.00949696, 0.0075684, 0.00557192, 0.0052629, 0.00665291,
        0.00960165, 0.00973791, 0.00920582, 0.0057934, 0.00709962,
        0.00623121, 0.00602675, 0.00842413, 0.00743056, 0.00662455,
        0.00550107, 0.00772382, 0.00673513, 0.00695548, 0.00655254,
        0.00693598, 0.0077793, 0.00507072, 0.00923823, 0.0096096,
        0.00775265, 0.00634011, 0.0099512, 0.00691597, 0.00846828,
        0.00844976, 0.00717155, 0.00599579, 0.0098329, 0.00531845,
        0.00742575, 0.00610365, 0.00646987, 0.00914264, 0.00683633,
        0.00541674, 0.00598155, 0.00930187, 0.00988514, 0.00633991,
        0.00837704, 0.00540599, 0.00861733, 0.00708218, 0.0095908,
        0.00655768, 0.00970733, 0.00751624, 0.00674446, 0.0082351,
        0.00624873, 0.00614882, 0.00598173, 0.0097995, 0.00746457,
        0.00875807, 0.00736996, 0.0079377, 0.00792069, 0.00989943,
        0.00834217, 0.00619885, 0.00507599, 0.00609341, 0.0072776,
        0.0069671, 0.00906163, 0.00892778, 0.00544548, 0.00976005,
        0.00763728, 0.00798202, 0.00702528, 0.0082475, 0.00935663,
        0.00836968, 0.00985049, 0.00850561, 0.0091086, 0.0052252,
        0.00836349, 0.00827376, 0.00550873, 0.00921194, 0.00807086,
        0.00549164, 0.00797234, 0.00739208, 0.00616647, 0.00509878,
        0.00682784, 0.00809926, 0.0066464, 0.00653627, 0.00875561,
        0.00879312, 0.00859383, 0.00550591, 0.00758083, 0.00778899,
        0.00872402, 0.00951589, 0.00684519, 0.00714332, 0.00866384,
        0.00831318, 0.00778935, 0.0067507, 0.00597676, 0.00591904,
        0.00540792, 0.005406, 0.00922899, 0.00691836, 0.0053037,
        0.00948213, 0.00611635, 0.00634062, 0.00597249, 0.00983751,
        0.0055627, 0.00861082, 0.00966044, 0.00834001, 0.00929363,
        0.00621224, 0.00836964, 0.00850436, 0.00729166, 0.00935273,
        0.00847193, 0.00947439, 0.00876602, 0.00760145, 0.00749344,
        0.00726864, 0.00510823, 0.00767571, 0.00711487, 0.00578767,
        0.00559535, 0.00724676, 0.00519957, 0.0099329, 0.0068906,
        0.00691055, 0.00525563, 0.00713336, 0.00507873, 0.00515047,
        0.0066955, 0.00910484, 0.00729411, 0.0050742, 0.0058161,
        0.00869961, 0.00869147, 0.00877261, 0.00675835, 0.00676138,
        0.00901038, 0.00699069, 0.00863596, 0.00790562, 0.00682171,
        0.00540003, 0.00558063, 0.00944779, 0.0072617, 0.00997002,
        0.00681948, 0.00624977, 0.0067527, 0.00671543, 0.00818678,
        0.00506369, 0.00881634, 0.00708207, 0.0071612, 0.00740558,
        0.00724606, 0.00748735, 0.00672952, 0.00726673, 0.00702326,
        0.00759121, 0.00811635, 0.0062052, 0.00754219, 0.00797311,
        0.00508474, 0.00760247, 0.00619647, 0.00702269, 0.00913265,
        0.00663118, 0.00741608, 0.00512371, 0.00654375, 0.00819861,
        0.00657581, 0.00602899, 0.00645328, 0.00977189, 0.00543401,
        0.00731679, 0.00529193, 0.00769329, 0.00573018, 0.00817042,
        0.00632199, 0.00845458, 0.00673573, 0.00502084, 0.00647447])
    period = 2.0
    transit_time = 0.5
    duration = 0.16
    depth = 0.2
    m = np.abs((t-transit_time+0.5*period) % period-0.5*period) < 0.5*duration
    y[m] = 1.0 - depth
    randn_arr = np.array([
        -1.00326528e-02, -8.45644428e-01,  9.11460610e-01, -1.37449688e+00,
        -5.47065645e-01, -7.55266106e-05, -1.21166803e-01, -2.00858547e+00,
        -9.20646543e-01,  1.68234342e-01, -1.31989156e+00,  1.26642930e+00,
        4.95180889e-01, -5.14240391e-01, -2.20292465e-01,  1.86156412e+00,
        9.35988451e-01,  3.80219145e-01, -1.41551877e+00,  1.62961132e+00,
        1.05240107e+00, -1.48405388e-01, -5.49698069e-01, -1.87903939e-01,
        -1.20193668e+00, -4.70785558e-01,  7.63160514e-01, -1.80762128e+00,
        -3.14074374e-01,  1.13755973e-01,  1.03568037e-01, -1.17893695e+00,
        -1.18215289e+00,  1.08916538e+00, -1.22452909e+00,  1.00865096e+00,
        -4.82365315e-01,  1.07979635e+00, -4.21078505e-01, -1.16647132e+00,
        8.56554856e-01, -1.73912222e-02,  1.44857659e+00,  8.92200085e-01,
        -2.29426629e-01, -4.49667602e-01,  2.33723433e-02,  1.90210018e-01,
        -8.81748527e-01,  8.41939573e-01, -3.97363492e-01, -4.23027745e-01,
        -5.40688337e-01,  2.31017267e-01, -6.92052602e-01,  1.34970110e-01,
        2.76660307e+00, -5.36094601e-02, -4.34004738e-01, -1.66768923e+00,
        5.02219248e-02, -1.10923094e+00, -3.75558119e-01,  1.51607594e-01,
        -1.73098945e+00,  1.57462752e-01,  3.04515175e-01, -1.29710002e+00,
        -3.92309192e-01, -1.83066636e+00,  1.57550094e+00,  3.30563277e-01,
        -1.79588501e-01, -1.63435831e-01,  1.13144361e+00, -9.41655519e-02,
        3.30816771e-01,  1.51862956e+00, -3.46167148e-01, -1.09263532e+00,
        -8.24500575e-01,  1.42866383e+00,  9.14283085e-02, -5.02331288e-01,
        9.73644380e-01,  9.97957386e-01, -4.75647768e-01, -9.71936837e-01,
        -1.57052860e+00, -1.79388892e+00, -2.64986452e-01, -8.93195947e-01,
        1.85847441e+00,  5.85377547e-02, -1.94214954e+00,  1.41872928e+00,
        1.61710309e-01,  7.04979480e-01,  6.82034777e-01,  2.96556567e-01,
        5.23342630e-01,  2.38760672e-01, -1.10638591e+00,  3.66732198e-01,
        1.02390550e+00, -2.10056413e-01,  5.51302218e-01,  4.19589145e-01,
        1.81565206e+00, -2.52750301e-01, -2.92004163e-01, -1.16931740e-01,
        -1.02391075e-01, -2.27261771e+00, -6.42609841e-01,  2.99885067e-01,
        -8.25651467e-03, -7.99339154e-01, -6.64779252e-01, -3.55613128e-01,
        -8.01571781e-01, -5.13050610e-01, -5.39390119e-01,  8.95370847e-01,
        1.01639127e+00,  9.33585094e-01,  4.26701799e-01, -7.08322484e-01,
        9.59830450e-01, -3.14250587e-01,  2.30522083e-02,  1.33822053e+00,
        8.39928561e-02,  2.47284030e-01, -1.41277949e+00,  4.87009294e-01,
        -9.80006647e-01,  1.01193966e+00, -1.84599177e-01, -2.23616884e+00,
        -3.58020103e-01, -2.28034538e-01,  4.85475226e-01,  6.70512391e-01,
        -3.27764245e-01,  1.01286819e+00, -3.16705533e+00, -7.13988998e-01,
        -1.11236427e+00, -1.25418351e+00,  9.59706371e-01,  8.29170399e-01,
        -7.75770020e-01,  1.17805700e+00,  1.01466892e-01, -4.21684101e-01,
        -6.92922796e-01, -7.78271726e-01,  4.72774857e-01,  6.50154901e-01,
        2.38501212e-01, -2.05021768e+00,  2.96358656e-01,  5.65396564e-01,
        -6.69205605e-01,  4.32505429e-02, -1.86388430e+00, -1.22996906e+00,
        -3.24235348e-01, -3.09751144e-01,  3.51679372e-01, -1.18692539e+00,
        -3.41206065e-01, -4.89779780e-01,  5.28010474e-01,  1.42104277e+00,
        1.72092032e+00, -1.56844005e+00, -4.80141918e-02, -1.11252931e+00,
        -6.47449515e-02,  4.22919280e-01,  8.14908987e-02, -4.90116988e-02,
        1.48303917e+00,  7.20989392e-01, -2.72654462e-01,  2.42113609e-02,
        8.70897807e-01,  6.09790506e-01, -4.25076104e-01, -1.77524284e+00,
        -1.18465749e+00,  1.45979225e-01, -1.78652685e+00, -1.52394498e-01,
        -4.53569176e-01,  9.99252803e-01, -1.31804382e+00, -1.93176898e+00,
        -4.19640742e-01,  6.34763132e-01,  1.06991860e+00, -9.09327017e-01,
        4.70263748e-01, -1.11143045e+00, -7.48827466e-01,  5.67594726e-01,
        7.18150543e-01, -9.99380749e-01,  4.74898323e-01, -1.86849981e+00,
        -2.02658907e-01, -1.13424803e+00, -8.07699340e-01, -1.27607735e+00,
        5.53626395e-01,  5.53874470e-01, -6.91200445e-01,  3.75582306e-01,
        2.61272553e-01, -1.28451754e-01,  2.15817020e+00, -8.40878617e-01,
        1.43050907e-02, -3.82387029e-01, -3.71780015e-01,  1.59412004e-01,
        -2.94395700e-01, -8.60426760e-01,  1.24227498e-01,  1.18233165e+00,
        9.42766380e-01,  2.03044488e-01, -7.35396814e-01,  1.86429600e-01,
        1.08464302e+00,  1.19118926e+00,  3.59687060e-01, -3.64357200e-01,
        -2.02752749e-01,  7.72045927e-01,  6.86346215e-01, -1.75769961e+00,
        6.58617565e-01,  7.11288340e-01, -8.87191425e-01, -7.64981116e-01,
        -7.57164098e-01, -6.80262803e-01, -1.41674959e+00,  3.13091930e-01,
        -7.85719399e-01, -7.03838361e-02, -4.97568783e-01,  2.55177521e-01,
        -1.01061704e+00,  2.45265375e-01,  3.89781016e-01,  8.27594585e-01,
        1.96776909e+00, -2.09210177e+00,  3.20314334e-01, -7.09162842e-01,
        -1.92505867e+00,  8.41630623e-01,  1.33219988e+00, -3.91627710e-01,
        2.10916296e-01, -6.40767402e-02,  4.34197668e-01,  8.80535749e-01,
        3.44937336e-01,  3.45769929e-01,  1.25973654e+00, -1.64662222e-01,
        9.23064571e-01, -8.22000422e-01,  1.60708495e+00,  7.37825392e-01,
        -4.03759534e-01, -2.11454815e+00, -3.10717131e-04, -1.18180941e+00,
        2.99634603e-01,  1.45116882e+00,  1.60059793e-01, -1.78012614e-01,
        3.42205404e-01,  2.85650196e-01, -2.36286411e+00,  2.40936864e-01,
        6.20277356e-01, -2.59341634e-01,  9.78559078e-01, -1.27674575e-01,
        7.66998762e-01,  2.27310511e+00, -9.63911290e-02, -1.94213217e+00,
        -3.36591724e-01, -1.72589000e+00,  6.11237826e-01,  1.30935097e+00,
        6.95879662e-01,  3.20308213e-01, -6.44925458e-01,  1.57564975e+00,
        7.53276212e-01,  2.84469557e-01,  2.04860319e-01,  1.11627359e-01,
        4.52216424e-01, -6.13327179e-01,  1.52524993e+00,  1.52339753e-01,
        6.00054450e-01, -4.33567278e-01,  3.74918534e-01, -2.28175243e+00,
        -1.11829888e+00, -3.14131532e-02, -1.32247311e+00,  2.43941406e+00,
        -1.66808131e+00,  3.45900749e-01,  1.65577315e+00,  4.81287059e-01,
        -3.10227553e-01, -5.52144084e-01,  6.73255489e-01, -8.00270681e-01,
        -1.19486110e-01,  6.91198606e-01, -3.07879027e-01,  8.75100102e-02,
        -3.04086293e-01, -9.69797604e-01,  1.18915048e+00,  1.39306624e+00,
        -3.16699954e-01, -2.65576159e-01, -1.77899339e-01,  5.38803274e-01,
        -9.05300265e-01, -8.85253056e-02,  2.62959055e-01,  6.42042149e-01,
        -2.78083727e+00,  4.03403210e-01,  3.45846762e-01,  1.00772824e+00,
        -5.26264015e-01, -5.18353205e-01,  1.20251659e+00, -1.56315671e+00,
        1.62909029e+00,  2.55589446e+00,  4.77451685e-01,  8.14098474e-01,
        -1.48958171e+00, -6.94559787e-01,  1.05786255e+00,  3.61815347e-01,
        -1.81427463e-01,  2.32869132e-01,  5.06976484e-01, -2.93095701e-01,
        -2.89459450e-02, -3.63073748e-02, -1.05227898e+00,  3.23594628e-01,
        1.80358591e+00,  1.73196213e+00, -1.47639930e+00,  5.70631220e-01,
        6.75503781e-01, -4.10510463e-01, -9.64200035e-01, -1.32081431e+00,
        -4.44703779e-01,  3.50009137e-01, -1.58058176e-01, -6.10933088e-01,
        -1.24915663e+00,  3.50716258e-01,  1.06654245e+00, -9.26921972e-01,
        4.48428964e-01, -1.87947524e+00, -6.57466109e-01,  7.29604120e-01,
        -1.11776721e+00, -6.04436725e-01,  1.41796683e+00, -7.32843980e-01,
        -8.53944819e-01,  5.75848362e-01,  1.95473356e+00, -2.39669947e-01,
        7.68735860e-01,  1.34576918e+00,  3.25552163e-01, -2.69917901e-01,
        -8.76326739e-01, -1.42521096e+00,  1.11170175e+00,  1.80957146e-01,
        1.33280094e+00,  9.88925316e-01, -6.16970520e-01, -1.18688670e+00,
        4.12669583e-01, -6.32506884e-01,  3.76689141e-01, -7.31151938e-01,
        -8.61225253e-01, -1.40990810e-01,  9.34100620e-01,  3.06539895e-01,
        1.17837515e+00, -1.23356170e+00, -1.05707714e+00, -8.91636992e-02,
        2.16570138e+00,  6.74286114e-01, -1.06661274e+00, -7.61404530e-02,
        2.20714791e-01, -5.68685746e-01,  6.13274991e-01, -1.56446138e-01,
        -2.99330718e-01,  1.26025679e+00, -1.70966090e+00, -9.61805342e-01,
        -8.17308981e-01, -8.47681070e-01, -7.28753045e-01,  4.88475958e-01,
        1.09653283e+00,  9.16041261e-01, -1.01956213e+00, -1.07417899e-01,
        4.52265213e-01,  2.40002952e-01,  1.30574740e+00, -6.75334236e-01,
        1.56319421e-01, -3.93230715e-01,  2.51075019e-01, -1.07889691e+00,
        -9.28937721e-01, -7.30110860e-01, -5.63669311e-01,  1.54792327e+00,
        1.17540191e+00, -2.12649671e-01,  1.72933294e-01, -1.59443602e+00,
        -1.79292347e-01,  1.59614713e-01,  1.14568421e+00,  3.26804720e-01,
        4.32890059e-01,  2.97762890e-01,  2.69001190e-01, -1.39675918e+00,
        -4.16757668e-01,  1.43488680e+00,  8.23896443e-01,  4.94234499e-01,
        6.67153092e-02,  6.59441396e-01, -9.44889409e-01, -1.58005956e+00,
        -3.82086552e-01,  5.37923058e-01,  1.07829882e-01,  1.01395868e+00,
        3.51450517e-01,  4.48421962e-02,  1.32748495e+00,  1.13237578e+00,
        -9.80913012e-02, -1.10304986e+00, -9.07361492e-01, -1.61451138e-01,
        -3.66811384e-01,  1.65776233e+00, -1.68013415e+00, -6.42577869e-02,
        -1.06622649e+00,  1.16801869e-01,  3.82264833e-01, -4.04896974e-01,
        5.30481414e-01, -1.98626941e-01, -1.79395613e-01, -4.17888725e-01])
    y += dy * randn_arr
    return t, y, dy, dict(period=period, transit_time=transit_time,
                          duration=duration, depth=depth)


def test_32bit_bug():
    rand = np.random.default_rng(42)
    t = rand.uniform(0, 10, 500)
    y = np.ones_like(t)
    y[np.abs((t + 1.0) % 2.0-1) < 0.08] = 1.0 - 0.1
    y += 0.01 * rand.standard_normal(len(t))

    model = BoxLeastSquares(t, y)
    results = model.autopower(0.16)
    assert_allclose(results.period[np.argmax(results.power)],
                    2.000412388152837)
    periods = np.linspace(1.9, 2.1, 5)
    results = model.power(periods, 0.16)
    assert_allclose(
        results.power,
        [0.01723948, 0.0643028, 0.1338783, 0.09428816, 0.03577543], rtol=1.1e-7)


@pytest.mark.parametrize("objective", ["likelihood", "snr"])
def test_correct_model(data, objective):
    t, y, dy, params = data
    model = BoxLeastSquares(t, y, dy)
    periods = np.exp(np.linspace(np.log(params["period"]) - 0.1,
                                 np.log(params["period"]) + 0.1, 1000))
    results = model.power(periods, params["duration"], objective=objective)
    ind = np.argmax(results.power)
    for k, v in params.items():
        assert_allclose(results[k][ind], v, atol=0.01)
    chi = (results.depth[ind]-params["depth"]) / results.depth_err[ind]
    assert np.abs(chi) < 1


@pytest.mark.parametrize("objective", ["likelihood", "snr"])
@pytest.mark.parametrize("offset", [False, True])
def test_fast_method(data, objective, offset):
    t, y, dy, params = data
    if offset:
        t = t - params["transit_time"] + params["period"]
    model = BoxLeastSquares(t, y, dy)
    periods = np.exp(np.linspace(np.log(params["period"]) - 1,
                                 np.log(params["period"]) + 1, 10))
    durations = params["duration"]
    results = model.power(periods, durations, objective=objective)

    assert_allclose_blsresults(results, model.power(periods, durations,
                                                    method="slow",
                                                    objective=objective))


def test_input_units(data):
    t, y, dy, params = data

    t_unit = u.day
    y_unit = u.mag

    with pytest.raises(u.UnitConversionError):
        BoxLeastSquares(t * t_unit, y * y_unit, dy * u.one)
    with pytest.raises(u.UnitConversionError):
        BoxLeastSquares(t * t_unit, y * u.one, dy * y_unit)
    with pytest.raises(u.UnitConversionError):
        BoxLeastSquares(t * t_unit, y, dy * y_unit)
    model = BoxLeastSquares(t*t_unit, y * u.one, dy)
    assert model.dy.unit == model.y.unit
    model = BoxLeastSquares(t*t_unit, y * y_unit, dy)
    assert model.dy.unit == model.y.unit
    model = BoxLeastSquares(t*t_unit, y*y_unit)
    assert model.dy is None


def test_period_units(data):
    t, y, dy, params = data
    t_unit = u.day
    y_unit = u.mag
    model = BoxLeastSquares(t * t_unit, y * y_unit, dy)

    p = model.autoperiod(params["duration"])
    assert p.unit == t_unit
    p = model.autoperiod(params["duration"] * 24 * u.hour)
    assert p.unit == t_unit
    with pytest.raises(u.UnitConversionError):
        model.autoperiod(params["duration"] * u.mag)

    p = model.autoperiod(params["duration"], minimum_period=0.5)
    assert p.unit == t_unit
    with pytest.raises(u.UnitConversionError):
        p = model.autoperiod(params["duration"], minimum_period=0.5*u.mag)

    p = model.autoperiod(params["duration"], maximum_period=0.5)
    assert p.unit == t_unit
    with pytest.raises(u.UnitConversionError):
        p = model.autoperiod(params["duration"], maximum_period=0.5*u.mag)

    p = model.autoperiod(params["duration"], minimum_period=0.5,
                         maximum_period=1.5)
    p2 = model.autoperiod(params["duration"], maximum_period=0.5,
                          minimum_period=1.5)
    assert_quantity_allclose(p, p2)


@pytest.mark.parametrize("method", ["fast", "slow"])
@pytest.mark.parametrize("with_err", [True, False])
@pytest.mark.parametrize("t_unit", [None, u.day])
@pytest.mark.parametrize("y_unit", [None, u.mag])
@pytest.mark.parametrize("objective", ["likelihood", "snr"])
def test_results_units(data, method, with_err, t_unit, y_unit, objective):
    t, y, dy, params = data

    periods = np.linspace(params["period"]-1.0, params["period"]+1.0, 3)

    if t_unit is not None:
        t = t * t_unit
    if y_unit is not None:
        y = y * y_unit
        dy = dy * y_unit
    if not with_err:
        dy = None

    model = BoxLeastSquares(t, y, dy)
    results = model.power(periods, params["duration"], method=method,
                          objective=objective)

    if t_unit is None:
        assert not has_units(results.period)
        assert not has_units(results.duration)
        assert not has_units(results.transit_time)
    else:
        assert results.period.unit == t_unit
        assert results.duration.unit == t_unit
        assert results.transit_time.unit == t_unit

    if y_unit is None:
        assert not has_units(results.power)
        assert not has_units(results.depth)
        assert not has_units(results.depth_err)
        assert not has_units(results.depth_snr)
        assert not has_units(results.log_likelihood)
    else:
        assert results.depth.unit == y_unit
        assert results.depth_err.unit == y_unit
        assert results.depth_snr.unit == u.one

        if dy is None:
            assert results.log_likelihood.unit == y_unit * y_unit
            if objective == "snr":
                assert results.power.unit == u.one
            else:
                assert results.power.unit == y_unit * y_unit
        else:
            assert results.log_likelihood.unit == u.one
            assert results.power.unit == u.one


def test_autopower(data):
    t, y, dy, params = data
    duration = params["duration"] + np.linspace(-0.1, 0.1, 3)

    model = BoxLeastSquares(t, y, dy)
    period = model.autoperiod(duration)
    results1 = model.power(period, duration)
    results2 = model.autopower(duration)

    assert_allclose_blsresults(results1, results2)


@pytest.mark.parametrize("with_units", [True, False])
def test_model(data, with_units):
    t, y, dy, params = data

    # Compute the model using linear regression
    A = np.zeros((len(t), 2))
    p = params["period"]
    dt = np.abs((t-params["transit_time"]+0.5*p) % p-0.5*p)
    m_in = dt < 0.5*params["duration"]
    A[~m_in, 0] = 1.0
    A[m_in, 1] = 1.0
    w = np.linalg.solve(np.dot(A.T, A / dy[:, None]**2),
                        np.dot(A.T, y / dy**2))
    model_true = np.dot(A, w)

    if with_units:
        t = t * u.day
        y = y * u.mag
        dy = dy * u.mag
        model_true = model_true * u.mag

    # Compute the model using the periodogram
    pgram = BoxLeastSquares(t, y, dy)
    model = pgram.model(t, p, params["duration"], params["transit_time"])

    # Make sure that the transit mask is consistent with the model
    transit_mask = pgram.transit_mask(t, p, params["duration"],
                                      params["transit_time"])
    transit_mask0 = (model - model.max()) < 0.0
    assert_allclose(transit_mask, transit_mask0)

    assert_quantity_allclose(model, model_true)


@pytest.mark.parametrize("shape", [(1,), (2,), (3,), (2, 3)])
def test_shapes(data, shape):
    t, y, dy, params = data
    duration = params["duration"]
    model = BoxLeastSquares(t, y, dy)

    period = np.empty(shape)
    period.flat = np.linspace(params["period"]-1, params["period"]+1,
                              period.size)
    if len(period.shape) > 1:
        with pytest.raises(ValueError):
            results = model.power(period, duration)
    else:
        results = model.power(period, duration)
        for k, v in results.items():
            if k == "objective":
                continue
            assert v.shape == shape


@pytest.mark.parametrize("with_units", [True, False])
@pytest.mark.parametrize("with_err", [True, False])
def test_compute_stats(data, with_units, with_err):
    t, y, dy, params = data

    y_unit = 1
    if with_units:
        y_unit = u.mag
        t = t * u.day
        y = y * u.mag
        dy = dy * u.mag
        params["period"] = params["period"] * u.day
        params["duration"] = params["duration"] * u.day
        params["transit_time"] = params["transit_time"] * u.day
        params["depth"] = params["depth"] * u.mag
    if not with_err:
        dy = None

    model = BoxLeastSquares(t, y, dy)
    results = model.power(params["period"], params["duration"],
                          oversample=1000)
    stats = model.compute_stats(params["period"], params["duration"],
                                params["transit_time"])

    # Test the calculated transit times
    tt = params["period"] * np.arange(int(t.max() / params["period"]) + 1)
    tt += params["transit_time"]
    assert_quantity_allclose(tt, stats["transit_times"])

    # Test that the other parameters are consistent with the periodogram
    assert_allclose(stats["per_transit_count"], [9, 7, 7, 7, 8])
    assert_quantity_allclose(np.sum(stats["per_transit_log_likelihood"]),
                             results["log_likelihood"])
    assert_quantity_allclose(stats["depth"][0], results["depth"])

    # Check the half period result
    results_half = model.power(0.5*params["period"], params["duration"],
                               oversample=1000)
    assert_quantity_allclose(stats["depth_half"][0], results_half["depth"])

    # Skip the uncertainty tests when the input errors are None
    if not with_err:
        assert_quantity_allclose(stats["harmonic_amplitude"],
                                 0.029945029964964204 * y_unit)
        assert_quantity_allclose(stats["harmonic_delta_log_likelihood"],
                                 -0.5875918155223113 * y_unit * y_unit)
        return

    assert_quantity_allclose(stats["harmonic_amplitude"],
                             0.033027988742275853 * y_unit)
    assert_quantity_allclose(stats["harmonic_delta_log_likelihood"],
                             -12407.505922833765)

    assert_quantity_allclose(stats["depth"][1], results["depth_err"])
    assert_quantity_allclose(stats["depth_half"][1], results_half["depth_err"])
    for f, k in zip((1.0, 1.0, 1.0, 0.0),
                    ("depth", "depth_even", "depth_odd", "depth_phased")):
        res = np.abs((stats[k][0]-f*params["depth"]) / stats[k][1])
        assert res < 1, f'f={f}, k={k}, res={res}'


def test_negative_times(data):
    t, y, dy, params = data
    mu = np.mean(t)
    duration = params["duration"] + np.linspace(-0.1, 0.1, 3)

    model1 = BoxLeastSquares(t, y, dy)
    results1 = model1.autopower(duration)

    # Compute the periodogram with offset (negative) times
    model2 = BoxLeastSquares(t - mu, y, dy)
    results2 = model2.autopower(duration)

    # Shift the transit times back into the unshifted coordinates
    results2.transit_time = (results2.transit_time + mu) % results2.period

    assert_allclose_blsresults(results1, results2)


@pytest.mark.parametrize('timedelta', [False, True])
def test_absolute_times(data, timedelta):

    # Make sure that we handle absolute times correctly. We also check that
    # TimeDelta works properly when timedelta is True.

    # The example data uses relative times
    t, y, dy, params = data

    # Add units
    t = t * u.day
    y = y * u.mag
    dy = dy * u.mag

    # We now construct a set of absolute times but keeping the rest the same.
    start = Time('2019-05-04T12:34:56')
    trel = TimeDelta(t) if timedelta else t
    t = trel + start

    # and we set up two instances of BoxLeastSquares, one with absolute and one
    # with relative times.
    bls1 = BoxLeastSquares(t, y, dy)
    bls2 = BoxLeastSquares(trel, y, dy)

    results1 = bls1.autopower(0.16 * u.day)
    results2 = bls2.autopower(0.16 * u.day)

    # All the results should match except transit time which should be
    # absolute instead of relative in the first case.

    for key in results1:
        if key == 'transit_time':
            assert_quantity_allclose((results1[key] - start).to(u.day), results2[key])
        elif key == 'objective':
            assert results1[key] == results2[key]
        else:
            assert_allclose(results1[key], results2[key])

    # Check that model evaluation works fine

    model1 = bls1.model(t, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    model2 = bls2.model(trel, 0.2 * u.day, 0.05 * u.day, TimeDelta(1 * u.day))
    assert_quantity_allclose(model1, model2)

    # Check model validation

    with pytest.raises(TypeError) as exc:
        bls1.model(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('transit_time was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls1.model(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('t_model was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls2.model(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('transit_time was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    with pytest.raises(TypeError) as exc:
        bls2.model(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('t_model was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    # Check compute_stats

    stats1 = bls1.compute_stats(0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    stats2 = bls2.compute_stats(0.2 * u.day, 0.05 * u.day, 1 * u.day)

    for key in stats1:
        if key == 'transit_times':
            assert_quantity_allclose((stats1[key] - start).to(u.day), stats2[key], atol=1e-10 * u.day)  # noqa: E501
        elif key.startswith('depth'):
            for value1, value2 in zip(stats1[key], stats2[key]):
                assert_quantity_allclose(value1, value2)
        else:
            assert_allclose(stats1[key], stats2[key])

    # Check compute_stats validation

    with pytest.raises(TypeError) as exc:
        bls1.compute_stats(0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('transit_time was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls2.compute_stats(0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('transit_time was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    # Check transit_mask

    mask1 = bls1.transit_mask(t, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    mask2 = bls2.transit_mask(trel, 0.2 * u.day, 0.05 * u.day, 1 * u.day)

    assert_equal(mask1, mask2)

    # Check transit_mask validation

    with pytest.raises(TypeError) as exc:
        bls1.transit_mask(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('transit_time was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls1.transit_mask(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('t was provided as a relative time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with absolute times.')

    with pytest.raises(TypeError) as exc:
        bls2.transit_mask(trel, 0.2 * u.day, 0.05 * u.day, Time('2019-06-04T12:34:56'))
    assert exc.value.args[0] == ('transit_time was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')

    with pytest.raises(TypeError) as exc:
        bls2.transit_mask(t, 0.2 * u.day, 0.05 * u.day, 1 * u.day)
    assert exc.value.args[0] == ('t was provided as an absolute time '
                                 'but the BoxLeastSquares class was initialized '
                                 'with relative times.')


def test_transit_time_in_range(data):
    t, y, dy, params = data

    t_ref = 10230.0
    t2 = t + t_ref
    bls1 = BoxLeastSquares(t, y, dy)
    bls2 = BoxLeastSquares(t2, y, dy)

    results1 = bls1.autopower(0.16)
    results2 = bls2.autopower(0.16)

    assert np.allclose(results1.transit_time, results2.transit_time - t_ref)
    assert np.all(results1.transit_time >= t.min())
    assert np.all(results1.transit_time <= t.max())
    assert np.all(results2.transit_time >= t2.min())
    assert np.all(results2.transit_time <= t2.max())
