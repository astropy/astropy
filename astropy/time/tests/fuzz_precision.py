#!/usr/bin/python
import atheris
import sys


with atheris.instrument_imports():
    import astropy.time.tests.test_precision as tp


def TestOneInput(data):
    fdp = atheris.FuzzedDataProvider(data)
    byte_size = fdp.ConsumeIntInRange(0, 1000)
    tp.test_abs_jd2_always_less_than_half_on_construction.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_round_to_even.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_two_sum.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_two_sum_symmetric.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_two_sum_size.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_day_frac_harmless.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_day_frac_exact.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_day_frac_idempotent.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_mjd_initialization_precise.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_day_frac_always_less_than_half.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_day_frac_round_to_even.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_resolution_never_decreases.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_resolution_never_decreases_utc.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_conversion_preserves_jd1_jd2_invariant.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_conversion_never_loses_precision.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_leap_stretch_mjd.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_jd_add_subtract_round_trip.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_time_argminmaxsort.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_timedelta_full_precision.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_timedelta_full_precision_arithmetic.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_timedelta_conversion.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_datetime_difference_agrees_with_timedelta.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_datetime_to_timedelta.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_datetime_timedelta_roundtrip.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_timedelta_datetime_roundtrip.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_timedelta_from_parts.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_datetime_timedelta_sum.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_sidereal_lat_independent.hypothesis.fuzz_one_input(
        fdp.ConsumeBytes(byte_size))
    tp.test_sidereal_lon_independent.hypothesis.fuzz_one_input(
       fdp.ConsumeBytes(byte_size))

tp.setup_module()
atheris.Setup(sys.argv, TestOneInput)
atheris.Fuzz()
