from distutils.core import Extension
from os.path import join
import sys

from astropy import setup_helpers


def get_build_options():
    return [
        ('use-system-expat',
         "Use the system expat library",
         True)
         ]


def get_extensions(build_type='release'):
    XML_DIR = 'astropy/utils/xml/src'

    source_files = [join(XML_DIR, "iterparse.c")]
    include_dirs = [XML_DIR]
    extra_link_args = []
    defines = []
    libraries = []
    library_dirs = []

    if setup_helpers.get_distutils_build_option('use_system_expat'):
        setup_helpers.pkg_config(
            ['expat'], ['expat'], include_dirs, library_dirs, libraries)
    else:
        EXPAT_DIR = 'cextern/expat/lib'
        source_files.extend([
            join(EXPAT_DIR, "xmlparse.c"),
            join(EXPAT_DIR, "xmlrole.c"),
            join(EXPAT_DIR, "xmltok.c"),
            join(EXPAT_DIR, "xmltok_impl.c")])
        include_dirs.append(EXPAT_DIR)
        if sys.platform.startswith('linux'):
            # This is to ensure we only export the Python entry point
            # symbols and the linker won't try to use the system expat in
            # place of ours.
            extra_link_args.append(
                '-Wl,--version-script={0}'.format(
                    join(XML_DIR, 'iterparse.map'))
                )
        defines = [("HAVE_EXPAT_CONFIG_H", 1)]
        if sys.byteorder == 'big':
            defines.append(('BYTEORDER', '4321'))
        else:
            defines.append(('BYTEORDER', '1234'))
        if sys.platform != 'win32':
            defines.append(('HAVE_UNISTD_H', None))

    return [Extension(
        "astropy.utils.xml._iterparser",
        source_files,
        define_macros=defines,
        include_dirs=include_dirs,
        library_dirs=library_dirs,
        libraries=libraries,
        extra_link_args=extra_link_args)]
