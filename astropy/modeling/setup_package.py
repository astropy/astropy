# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
from os.path import join
from collections import defaultdict

from setuptools import Extension

from extension_helpers import import_file

wcs_setup_package = import_file(join('astropy', 'wcs', 'setup_package.py'))


MODELING_ROOT = os.path.relpath(os.path.dirname(__file__))
MODELING_SRC = join(MODELING_ROOT, 'src')
SRC_FILES = [join(MODELING_SRC, 'projections.c.templ'),
             __file__]
GEN_FILES = [join(MODELING_SRC, 'projections.c')]


# This defines the set of projection functions that we want to wrap.
# The key is the projection name, and the value is the number of
# parameters.

# (These are in the order that the appear in the WCS coordinate
# systems paper).
projections = {
    'azp': 2,
    'szp': 3,
    'tan': 0,
    'stg': 0,
    'sin': 2,
    'arc': 0,
    'zea': 0,
    'air': 1,
    'cyp': 2,
    'cea': 1,
    'mer': 0,
    'sfl': 0,
    'par': 0,
    'mol': 0,
    'ait': 0,
    'cop': 2,
    'coe': 2,
    'cod': 2,
    'coo': 2,
    'bon': 1,
    'pco': 0,
    'tsc': 0,
    'csc': 0,
    'qsc': 0,
    'hpx': 2,
    'xph': 0,
}


def get_extensions():

    from jinja2 import Environment, FileSystemLoader

    # Prepare the jinja2 templating environment
    env = Environment(loader=FileSystemLoader(MODELING_SRC))

    c_in = env.get_template('projections.c.templ')
    c_out = c_in.render(projections=projections)

    with open(join(MODELING_SRC, 'projections.c'), 'w') as fd:
        fd.write(c_out)

    wcslib_files = [  # List of wcslib files to compile
        'prj.c',
        'wcserr.c',
        'wcsprintf.c',
        'wcsutil.c'
    ]

    wcslib_config_paths = [
        join(MODELING_SRC, 'wcsconfig.h')
    ]

    cfg = defaultdict(list)

    wcs_setup_package.get_wcslib_cfg(cfg, wcslib_files, wcslib_config_paths)

    cfg['include_dirs'].append(MODELING_SRC)

    astropy_files = [  # List of astropy.modeling files to compile
        'projections.c'
    ]
    cfg['sources'].extend(join(MODELING_SRC, x) for x in astropy_files)

    cfg['sources'] = [str(x) for x in cfg['sources']]
    cfg = dict((str(key), val) for key, val in cfg.items())

    return [Extension('astropy.modeling._projections', **cfg)]
