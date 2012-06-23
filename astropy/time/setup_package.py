import os
import numpy
from distutils.extension import Extension

TIMEROOT = os.path.relpath(os.path.dirname(__file__))


def get_extensions():
    time_ext = Extension(
    name="astropy.time.sofa_time",
    sources=[os.path.join(TIMEROOT, x) for x in ("sofa_time.pyx", "sofa.c")],
    include_dirs=[numpy.get_include()],
    language="c",)
    return [time_ext]
