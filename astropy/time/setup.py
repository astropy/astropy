# python f-setup.py build_ext --inplace
#   cython f.pyx -> f.cpp
#   g++ -c f.cpp -> f.o
#   g++ -c fc.cpp -> fc.o
#   link f.o fc.o -> f.so

# distutils uses the Makefile distutils.sysconfig.get_makefile_filename()
import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name="sofa_time",
    sources=["sofa_time.pyx", "sofa.c"],
    include_dirs = [numpy.get_include()],
    language="c",
    )]

setup(
    name = 'astrotime',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    )
