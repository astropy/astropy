# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Build a small C extension against the public ``astropy.wcs`` C API.

The extension links against the headers published by ``astropy.wcs.get_include()``,
exactly as a downstream package such as drizzlepac does.  It is used to check
that the members of the C API still in downstream use keep working when called
through the function-pointer table, and that the deprecated members emit a
``DeprecationWarning``.
"""

import importlib
import os
import shutil
import subprocess
import sys
import sysconfig
import warnings

import pytest

from astropy.wcs import WCS

HERE = os.path.dirname(__file__)

pytestmark = pytest.mark.skipif(
    not os.path.exists(os.path.join(sysconfig.get_path("include"), "Python.h")),
    reason="CPython development headers are required to build the test extension",
)

# Built at test time against the headers published by astropy.wcs.get_include(),
# with the same include directories a downstream package such as drizzlepac uses.
SETUP_PY = """\
import os
import sys

import numpy as np
from setuptools import Extension, setup

from astropy import wcs

if sys.platform == "win32":
    define_macros = [
        ("YY_NO_UNISTD_H", None),
        ("_CRT_SECURE_NO_WARNINGS", None),
        ("_NO_OLDNAMES", None),
        ("NO_OLDNAMES", None),
        ("__STDC__", None),
    ]
else:
    define_macros = []

setup(
    name="wcsapi_test",
    ext_modules=[
        Extension(
            "wcsapi_test",
            sources=["wcsapi_test.c"],
            include_dirs=[
                np.get_include(),
                os.path.join(wcs.get_include(), "astropy_wcs"),
                os.path.join(wcs.get_include(), "wcslib"),
            ],
            define_macros=define_macros,
        )
    ],
)
"""


@pytest.fixture(scope="module")
def wcsapi_test(tmp_path_factory):
    # Build the extension in a scratch directory so we don't leave artifacts in
    # the source tree.
    build_dir = tmp_path_factory.mktemp("wcsapi_test")
    shutil.copy(os.path.join(HERE, "wcsapi_test.c"), build_dir)
    (build_dir / "setup.py").write_text(SETUP_PY)

    proc = subprocess.run(
        [sys.executable, "setup.py", "build_ext", "--inplace"],
        cwd=build_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        message = (
            "building the wcsapi_test extension failed:\n"
            f"{proc.stdout}\n{proc.stderr}"
        )
        # On CI we want to know if this breaks; elsewhere a missing compiler or
        # development headers should not fail the suite.
        if os.environ.get("CI") or os.environ.get("CONTINUOUS_INTEGRATION"):
            raise AssertionError(message)
        pytest.skip(message)

    sys.path.insert(0, str(build_dir))
    try:
        yield importlib.import_module("wcsapi_test")
    finally:
        sys.path.remove(str(build_dir))
        sys.modules.pop("wcsapi_test", None)


def _make_wcs():
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, 20.0]
    w.wcs.crpix = [128.0, 128.0]
    w.wcs.cdelt = [-0.001, 0.001]
    w.wcs.set()
    return w


def test_get_c_version(wcsapi_test):
    # Slot 0: the version read back through the table is the one import checked.
    assert wcsapi_test.get_c_version() > 0


def test_error_message(wcsapi_test):
    # Slot 23: wcslib_get_error_message returns a string for a known code.
    assert wcsapi_test.get_error_message(1) != ""


def test_roundtrip(wcsapi_test):
    # Slots 1, 2, 20, 21: wcsprm_python2c, wcsprm_c2python, wcsp2s and wcss2p.
    w = _make_wcs()
    x, y = 200.0, 150.0
    rx, ry = wcsapi_test.roundtrip(w.wcs, x, y)
    assert rx == pytest.approx(x, abs=1e-6)
    assert ry == pytest.approx(y, abs=1e-6)


def test_print_wcsprm(wcsapi_test):
    # Slot 22: wcsprt returns 0 for a valid, set wcsprm.
    w = _make_wcs()
    assert wcsapi_test.print_wcsprm(w.wcs) == 0


def test_all_pix2world(wcsapi_test):
    # Slot 18: pipeline_all_pixel2world agrees with WCS.all_pix2world.
    w = _make_wcs()
    x, y, origin = 200.0, 150.0, 1
    cx, cy = wcsapi_test.all_pix2world(w, x, y, origin)
    (ex, ey),  = w.all_pix2world([[x, y]], origin)
    assert cx == pytest.approx(ex, abs=1e-9)
    assert cy == pytest.approx(ey, abs=1e-9)


def test_supported_slots_do_not_warn(wcsapi_test):
    # None of the kept members of the table warn.
    w = _make_wcs()
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        wcsapi_test.get_c_version()
        wcsapi_test.get_error_message(1)
        wcsapi_test.roundtrip(w.wcs, 200.0, 150.0)
        wcsapi_test.print_wcsprm(w.wcs)
        wcsapi_test.all_pix2world(w, 200.0, 150.0, 1)


def test_deprecated_slot_warns(wcsapi_test):
    # A deprecated member of the table (pipeline_clear) warns when used.
    with pytest.warns(DeprecationWarning, match="pipeline_clear"):
        wcsapi_test.call_deprecated()
