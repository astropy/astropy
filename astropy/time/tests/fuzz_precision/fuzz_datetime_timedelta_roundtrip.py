#!/usr/bin/python
import atheris
import sys


with atheris.instrument_imports():
    import astropy.time.tests.test_precision as tp

tp.setup_module()
atheris.Setup(
    sys.argv, tp.test_datetime_timedelta_roundtrip.hypothesis.fuzz_one_input)
atheris.Fuzz()
