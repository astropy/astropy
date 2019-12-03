# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from distutils.core import Extension
from os.path import join
import sys

from astropy_helpers import setup_helpers


def get_extensions(build_type='release'):
    XML_DIR = 'astropy/utils/xml/src'

    cfg = setup_helpers.DistutilsExtensionArgs({
        'sources': [join(XML_DIR, "iterparse.c")]
        })

    if int(os.environ.get('ASTROPY_USE_SYSTEM_EXPAT', 0)) or int(os.environ.get('ASTROPY_USE_SYSTEM_ALL', 0)):
        cfg.update(setup_helpers.pkg_config(['expat'], ['expat']))
    else:
        EXPAT_DIR = 'cextern/expat/lib'
        cfg['sources'].extend([
            join(EXPAT_DIR, fn) for fn in
            ["xmlparse.c", "xmlrole.c", "xmltok.c", "xmltok_impl.c",
             "loadlibrary.c"]])
        cfg['include_dirs'].extend([XML_DIR, EXPAT_DIR])
        if sys.platform.startswith('linux'):
            # This is to ensure we only export the Python entry point
            # symbols and the linker won't try to use the system expat in
            # place of ours.
            cfg['extra_link_args'].extend([
                '-Wl,--version-script={}'.format(
                    join(XML_DIR, 'iterparse.map'))
                ])
        cfg['define_macros'].append(("HAVE_EXPAT_CONFIG_H", 1))
        if sys.byteorder == 'big':
            cfg['define_macros'].append(('BYTEORDER', '4321'))
        else:
            cfg['define_macros'].append(('BYTEORDER', '1234'))
        if sys.platform != 'win32':
            cfg['define_macros'].append(('HAVE_UNISTD_H', None))

    return [Extension("astropy.utils.xml._iterparser", **cfg)]
