# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file needs to be included here to make sure commands such
# as ``python setup.py test ... -t docs/...`` works, since this
# will ignore the conftest.py file at the root of the repository
# and the one in astropy/conftest.py

import os
import tempfile

# Make sure we use temporary directories for the config and cache
# so that the tests are insensitive to local configuration.

os.environ['XDG_CONFIG_HOME'] = tempfile.mkdtemp('astropy_config')
os.environ['XDG_CACHE_HOME'] = tempfile.mkdtemp('astropy_cache')

os.mkdir(os.path.join(os.environ['XDG_CONFIG_HOME'], 'astropy'))
os.mkdir(os.path.join(os.environ['XDG_CACHE_HOME'], 'astropy'))

# Note that we don't need to change the environment variables back or remove
# them after testing, because they are only changed for the duration of the
# Python process, and this configuration only matters if running pytest
# directly, not from e.g. an IPython session.
