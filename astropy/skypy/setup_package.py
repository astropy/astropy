#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup
from setuptools.extension import Extension
import os
import glob

from .. import __path__ as astropy_path
# Following two lines needed for the alternative method described at
# the end.
# from .. import setup_helpers
# from ..version import release

# These we would in a skypy.cfg and there can be defaults for
# paths. Here the target header is ~/lib/skysub.h and library is
# ~/lib/libskycal.a.
USE_SIN_FUNC_LIB = False
SIN_FUNC_LIB_PATH = os.path.abspath(os.path.join(
        os.path.expanduser("~"), "lib"))
SIN_FUNC_INCLUDE_PATH = SIN_FUNC_LIB_PATH

# Find skycalc source directory shipped with astropy.
# cextern/skycalc.
cextern_dir = os.path.abspath(os.path.join(
        os.path.dirname(astropy_path[0]), "cextern"
        ))
src_dir = os.path.abspath(os.path.join(
        cextern_dir, "skycalc"))

this_dir = os.path.abspath(os.path.dirname(__file__))
# This package has two extensions: skypy.skyc and skypy.astrom. The
# first uses external C library, and the other is a stand alone Cython
# code. The former can either use the C source shipped with astropy, or
# directly link with library on user's machine.

# Extension skypy.skyc. Build the Cython file, .pyx or .c, along with
# the skycalc.c file shipped with astropy, and create a skyc.so
# file. Or Build the Cython file, .pyx or .c, and then create skyc.so
# that links with user's libskycalc.
src_files = []
depends = []
include_dirs = []
libraries = []
library_dirs = []

if USE_SIN_FUNC_LIB:
    # Don't use skycalc that comes with astropy.
    library_dirs.append(SIN_FUNC_LIB_PATH)
    include_dirs.append(SIN_FUNC_INCLUDE_PATH)
    libraries.append('skycalc')
else:
    # Compile sykc.c, or sky.c from skyc.pyx, along with the skycalc.c
    # file shipped with astropy.
    src_files = glob.glob(os.path.join(src_dir, "*.c"))
    depends = glob.glob(os.path.join(src_dir, "*.h"))
    depends.extend(glob.glob(os.path.join(this_dir, "skyc.pxd")))
    # This could be src_dir or cextern_dir based on cdef extern from "<x>".
    include_dirs = [src_dir]


# setup.py will call this function if ``HAVE_CYTHON and not release``,
# otherwise it will call the next function. setup.py should call these
# and not try to find Cython extensions by itself, for the same reason
# that setup.py doesn't try to find regular C extensions by itself.
def get_cython_pyx_extensions():
    # Run Cython and create skyc.c from skyc.pyx.
    src_files.append(os.path.join(this_dir, "skyc.pyx"))
    skyc_ext = Extension('astropy.skypy.skyc',
                         src_files,
                         depends=depends,
                         include_dirs=include_dirs,
                         library_dirs=library_dirs,
                         libraries=libraries)
    astrom_ext = Extension('astropy.skypy.astrom',
                           [os.path.join(this_dir, 'astrom.pyx')])

    return [skyc_ext, astrom_ext]


def get_cython_extensions():
    # Use existing skyc.c without running Cython.
    src_files.append(os.path.join(this_dir, "skyc.c"))
    skyc_ext = Extension('astropy.skypy.skyc',
                         src_files,
                         depends=depends,
                         include_dirs=include_dirs,
                         library_dirs=library_dirs,
                         libraries=libraries)
    astrom_ext = Extension('astropy.skypy.astrom',
                           [os.path.join(this_dir, 'astrom.c')])

    return [skyc_ext, astrom_ext]


# Alternative method, where developer checks if pyx files or c files are
#to be included and then returns appropriate values in
#get_extension(). In this scenario there is not need for a
#get_cython_extensions() or get_cython_pyx_extension() function.
#if setup_helpers.HAVE_CYTHON and not release:
#    # Run Cython and create skyc.c from skyc.pyx.
#    src_files.append(os.path.join("skyc.pyx"))
#else:
#    # Use existing skyc.c without running Cython.
#    src_files.append(os.path.join("skyc.c"))
#
#skyc_ext = Extension('astropy.skypy.skyc',
#                     src_files,
#                     depends=depends,
#                     include_dirs=include_dirs,
#                     library_dirs=library_dirs,
#                     libraries=libraries)
#
## A simple extension that doesn't use any external C code.
#if setup_helpers.HAVE_CYTHON and not release:
#    astrom_ext = Extension('astropy.skypy.astrom', ['astrom.pyx'])
#else:
#    astrom_ext = Extension('astropy.skypy.astrom', ['astrom.c'])
#
# # This is always called by setup.py.
# def get_extension():
#    # Add pure C extensions as well.
#    return [skyc_ext, astrom_ext]
