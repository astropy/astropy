# Licensed under a 3-clause BSD style license - see LICENSE.rst

from os.path import join
import sysconfig

from setuptools import Extension

USE_PY_LIMITED_API = not sysconfig.get_config_var("Py_GIL_DISABLED")


def get_extensions(build_type="release"):
    VO_DIR = "astropy/io/votable/src"

    kwargs = {}
    if USE_PY_LIMITED_API:
        kwargs["py_limited_api"] = True
        kwargs['define_macros'] = [("Py_LIMITED_API", "0x030B0000")]

    return [
        Extension(
            "astropy.io.votable.tablewriter",
            [join(VO_DIR, "tablewriter.c")],
            include_dirs=[VO_DIR],
            **kwargs
        )
    ]
