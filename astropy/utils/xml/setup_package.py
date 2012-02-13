from distutils.core import Extension
from os.path import join
import sys


def get_extensions(build_type='release'):
    EXPAT_DIR = 'cextern/expat/lib'
    XML_DIR = 'astropy/utils/xml/src'

    defines = [("HAVE_EXPAT_CONFIG_H", 1)]
    if sys.byteorder == 'big':
        defines.append(('BYTEORDER', '4321'))
    else:
        defines.append(('BYTEORDER', '1234'))
    if sys.platform != 'win32':
        defines.append(('HAVE_UNISTD_H', None))

    if sys.platform.startswith('linux'):
        # This is to ensure we only export the Python entry point
        # symbols and the linker won't try to use the system expat in
        # place of ours.
        extra_link_args = [
            '-Wl,--version-script={0}'.format(
                join(XML_DIR, 'iterparse.map'))
            ]
    else:
        extra_link_args = []

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
