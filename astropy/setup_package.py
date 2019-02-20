# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import glob


def get_package_data():

    # Find all files in data/ sub-directories since this is a standard location
    # for data files. We then need to adjust the paths to be relative to here
    # (otherwise glob will be evaluated relative to setup.py)
    data_files = glob.glob('**/data/**/*', recursive=True)
    data_files = [os.path.relpath(x, os.path.dirname(__file__)) for x in data_files]

    # Glob doesn't recognize hidden files
    data_files.append('utils/tests/data/.hidden_file.txt')

    return {'astropy': ['astropy.cfg', 'CITATION'] + data_files}
