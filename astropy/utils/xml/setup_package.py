from distutils.core import Extension
from os.path import join
import sys


def get_extensions(build_type='release'):
    EXPAT_DIR = 'cextern/expat/lib'
    XML_DIR = 'astropy/utils/xml/src'

    if sys.platform.startswith('linux'):
        # This is to ensure that _iterparser uses the libexpat functions it was
        # compiled with, and not the system libexpat (which can happen if the
        # latter is loaded first).
        extra_link_args = ['-Wl,-Bsymbolic']
    else:
        extra_link_args = []

    defines = [("HAVE_EXPAT_CONFIG_H", 1)]
    if sys.byteorder == 'big':
        defines.append(('BYTEORDER', '4321'))
    else:
        defines.append(('BYTEORDER', '1234'))
    if sys.platform != 'win32':
        defines.append(('HAVE_UNISTD_H', None))

    return [Extension(
        "astropy.utils.xml._iterparser",
        [join(XML_DIR, "iterparse.c"),
         join(EXPAT_DIR, "xmlparse.c"),
         join(EXPAT_DIR, "xmlrole.c"),
         join(EXPAT_DIR, "xmltok.c"),
         join(EXPAT_DIR, "xmltok_impl.c")],
        define_macros=defines,
        include_dirs=[XML_DIR, EXPAT_DIR],
        extra_link_args=extra_link_args)]
