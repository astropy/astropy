# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

CONTACT = "Michael Droettboom"
EMAIL = "mdroe@stsci.edu"

import io
from os.path import join
import os.path
import shutil
import sys

from distutils.core import Extension
from distutils.dep_util import newer_group


from astropy_helpers import setup_helpers
from astropy_helpers.distutils_helpers import get_distutils_build_option
from astropy.extern import six

WCSROOT = os.path.relpath(os.path.dirname(__file__))
WCSVERSION = "5.10"


def b(s):
    return s.encode('ascii')

if six.PY2:
    def string_escape(s):
        # string_escape has subtle differences with the escaping done in Python
        # 3 so correct for those too
        s = s.encode('string_escape')
        s = s.replace(r'\x00', r'\0')
        return s.replace(r"\'", "'")
else:
    def string_escape(s):
        s = s.decode('ascii').encode('ascii', 'backslashreplace')
        s = s.replace(b'\n', b'\\n')
        s = s.replace(b'\0', b'\\0')
        return s.decode('ascii')


def determine_64_bit_int():
    """
    The only configuration parameter needed at compile-time is how to
    specify a 64-bit signed integer.  Python's ctypes module can get us
    that information, but it is only available in Python 2.5 or later.
    If we can't be absolutely certain, we default to "long long int",
    which is correct on most platforms (x86, x86_64).  If we find
    platforms where this heuristic doesn't work, we may need to
    hardcode for them.
    """
    try:
        try:
            import ctypes
        except ImportError:
            raise ValueError()

        if ctypes.sizeof(ctypes.c_longlong) == 8:
            return "long long int"
        elif ctypes.sizeof(ctypes.c_long) == 8:
            return "long int"
        elif ctypes.sizeof(ctypes.c_int) == 8:
            return "int"
        else:
            raise ValueError()

    except ValueError:
        return "long long int"


def write_wcsconfig_h(paths):
    """
    Writes out the wcsconfig.h header with local configuration.
    """
    h_file = io.StringIO()
    h_file.write("""
    /* The bundled version has WCSLIB_VERSION */
    #define HAVE_WCSLIB_VERSION 1

    /* WCSLIB library version number. */
    #define WCSLIB_VERSION {0}

    /* 64-bit integer data type. */
    #define WCSLIB_INT64 {1}

    /* Windows needs some other defines to prevent inclusion of wcsset()
       which conflicts with wcslib's wcsset().  These need to be set
       on code that *uses* astropy.wcs, in addition to astropy.wcs itself.
       */
    #if defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__) || defined (__MINGW64__)

    #ifndef YY_NO_UNISTD_H
    #define YY_NO_UNISTD_H
    #endif

    #ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS
    #endif

    #ifndef _NO_OLDNAMES
    #define _NO_OLDNAMES
    #endif

    #ifndef NO_OLDNAMES
    #define NO_OLDNAMES
    #endif

    #ifndef __STDC__
    #define __STDC__ 1
    #endif

    #endif
    """.format(WCSVERSION, determine_64_bit_int()))
    content = h_file.getvalue().encode('ascii')
    for path in paths:
        setup_helpers.write_if_different(path, content)


######################################################################
# GENERATE DOCSTRINGS IN C


def generate_c_docstrings():
    from astropy.wcs import docstrings
    docstrings = docstrings.__dict__
    keys = [
        key for key, val in docstrings.items()
        if not key.startswith('__') and isinstance(val, six.string_types)]
    keys.sort()
    docs = {}
    for key in keys:
        docs[key] = docstrings[key].encode('utf8').lstrip() + b'\0'

    h_file = io.StringIO()
    h_file.write("""/*
DO NOT EDIT!

This file is autogenerated by astropy/wcs/setup_package.py.  To edit
its contents, edit astropy/wcs/docstrings.py
*/

#ifndef __DOCSTRINGS_H__
#define __DOCSTRINGS_H__

""")
    for key in keys:
        val = docs[key]
        h_file.write('extern char doc_{0}[{1}];\n'.format(key, len(val)))
    h_file.write("\n#endif\n\n")

    setup_helpers.write_if_different(
        join(WCSROOT, 'include', 'astropy_wcs', 'docstrings.h'),
        h_file.getvalue().encode('utf-8'))

    c_file = io.StringIO()
    c_file.write("""/*
DO NOT EDIT!

This file is autogenerated by astropy/wcs/setup_package.py.  To edit
its contents, edit astropy/wcs/docstrings.py

The weirdness here with strncpy is because some C compilers, notably
MSVC, do not support string literals greater than 256 characters.
*/

#include <string.h>
#include "astropy_wcs/docstrings.h"

""")
    for key in keys:
        val = docs[key]
        c_file.write('char doc_{0}[{1}] = {{\n'.format(key, len(val)))
        for i in range(0, len(val), 12):
            section = val[i:i+12]
            if six.PY2:
                section = [ord(x) for x in section]
            c_file.write('    ');
            c_file.write(''.join('0x{0:02x}, '.format(x) for x in section))
            c_file.write('\n')

        c_file.write("    };\n\n")

    setup_helpers.write_if_different(
        join(WCSROOT, 'src', 'docstrings.c'),
        c_file.getvalue().encode('utf-8'))


def get_wcslib_cfg(cfg, wcslib_files, include_paths):
    from astropy.version import debug

    cfg['include_dirs'].append('numpy')
    cfg['define_macros'].extend([
        ('ECHO', None),
        ('WCSTRIG_MACRO', None),
        ('ASTROPY_WCS_BUILD', None),
        ('_GNU_SOURCE', None)])

    if (not setup_helpers.use_system_library('wcslib') or
            sys.platform == 'win32'):
        write_wcsconfig_h(include_paths)

        wcslib_path = join("cextern", "wcslib")  # Path to wcslib
        wcslib_cpath = join(wcslib_path, "C")  # Path to wcslib source files
        cfg['sources'].extend(join(wcslib_cpath, x) for x in wcslib_files)
        cfg['include_dirs'].append(wcslib_cpath)
    else:
        wcsconfig_h_path = join(WCSROOT, 'include', 'wcsconfig.h')
        if os.path.exists(wcsconfig_h_path):
            os.unlink(wcsconfig_h_path)
        cfg.update(setup_helpers.pkg_config(['wcslib'], ['wcs']))

    if debug:
        cfg['define_macros'].append(('DEBUG', None))
        cfg['undef_macros'].append('NDEBUG')
        if (not sys.platform.startswith('sun') and
            not sys.platform == 'win32'):
            cfg['extra_compile_args'].extend(["-fno-inline", "-O0", "-g"])
    else:
        # Define ECHO as nothing to prevent spurious newlines from
        # printing within the libwcs parser
        cfg['define_macros'].append(('NDEBUG', None))
        cfg['undef_macros'].append('DEBUG')

    if sys.platform == 'win32':
        # These are written into wcsconfig.h, but that file is not
        # used by all parts of wcslib.
        cfg['define_macros'].extend([
            ('YY_NO_UNISTD_H', None),
            ('_CRT_SECURE_NO_WARNINGS', None),
            ('_NO_OLDNAMES', None),  # for mingw32
            ('NO_OLDNAMES', None),  # for mingw64
            ('__STDC__', None)  # for MSVC
        ])

    if sys.platform.startswith('linux'):
        cfg['define_macros'].append(('HAVE_SINCOS', None))

    # Squelch a few compilation warnings in WCSLIB
    if setup_helpers.get_compiler_option() in ('unix', 'mingw32'):
        if not get_distutils_build_option('debug'):
            cfg['extra_compile_args'].extend([
                '-Wno-strict-prototypes',
                '-Wno-unused-function',
                '-Wno-unused-value',
                '-Wno-uninitialized'])



def get_extensions():
    generate_c_docstrings()

    ######################################################################
    # DISTUTILS SETUP
    cfg = setup_helpers.DistutilsExtensionArgs()

    wcslib_files = [  # List of wcslib files to compile
        'flexed/wcsbth.c',
        'flexed/wcspih.c',
        'flexed/wcsulex.c',
        'flexed/wcsutrn.c',
        'cel.c',
        'dis.c',
        'lin.c',
        'log.c',
        'prj.c',
        'spc.c',
        'sph.c',
        'spx.c',
        'tab.c',
        'wcs.c',
        'wcserr.c',
        'wcsfix.c',
        'wcshdr.c',
        'wcsprintf.c',
        'wcsunits.c',
        'wcsutil.c'
    ]

    wcslib_config_paths = [
        join(WCSROOT, 'include', 'astropy_wcs', 'wcsconfig.h'),
        join(WCSROOT, 'include', 'wcsconfig.h')
    ]

    get_wcslib_cfg(cfg, wcslib_files, wcslib_config_paths)

    cfg['include_dirs'].append(join(WCSROOT, "include"))

    astropy_wcs_files = [  # List of astropy.wcs files to compile
        'distortion.c',
        'distortion_wrap.c',
        'docstrings.c',
        'pipeline.c',
        'pyutil.c',
        'astropy_wcs.c',
        'astropy_wcs_api.c',
        'sip.c',
        'sip_wrap.c',
        'str_list_proxy.c',
        'unit_list_proxy.c',
        'util.c',
        'wcslib_wrap.c',
        'wcslib_tabprm_wrap.c']
    cfg['sources'].extend(join(WCSROOT, 'src', x) for x in astropy_wcs_files)

    cfg['sources'] = [str(x) for x in cfg['sources']]
    cfg = dict((str(key), val) for key, val in six.iteritems(cfg))

    return [Extension(str('astropy.wcs._wcs'), **cfg)]


def get_package_data():
    # Installs the testing data files
    api_files = [
        'astropy_wcs.h',
        'astropy_wcs_api.h',
        'distortion.h',
        'isnan.h',
        'pipeline.h',
        'pyutil.h',
        'sip.h',
        'util.h',
        'wcsconfig.h',
        ]
    api_files = [join('include', 'astropy_wcs', x) for x in api_files]
    api_files.append(join('include', 'astropy_wcs_api.h'))

    wcslib_headers = [
        'cel.h',
        'lin.h',
        'prj.h',
        'spc.h',
        'spx.h',
        'tab.h',
        'wcs.h',
        'wcserr.h',
        'wcsmath.h',
        'wcsprintf.h',
        ]
    if not setup_helpers.use_system_library('wcslib'):
        for header in wcslib_headers:
            source = join('cextern', 'wcslib', 'C', header)
            dest = join('astropy', 'wcs', 'include', 'wcslib', header)
            if newer_group([source], dest, 'newer'):
                shutil.copy(source, dest)
            api_files.append(join('include', 'wcslib', header))

    return {
        str('astropy.wcs.tests'): ['data/*.hdr', 'data/*.fits',
                                   'data/*.txt', 'data/*.fits.gz',
                                   'maps/*.hdr', 'spectra/*.hdr',
                                   'extension/*.c'],
        str('astropy.wcs'): api_files,
    }


def get_external_libraries():
    return ['wcslib']


def requires_2to3():
    return False
