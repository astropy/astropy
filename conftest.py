# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file is the main file used when running tests with pytest directly,
# in particular if running e.g. ``pytest docs/``.

import os
import tempfile

import hypothesis

from astropy import __version__

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
except ImportError:
    PYTEST_HEADER_MODULES = {}
    TESTED_VERSIONS = {}


# This has to be in the root dir or it will not display in CI.
def pytest_configure(config):
    PYTEST_HEADER_MODULES['PyERFA'] = 'erfa'
    PYTEST_HEADER_MODULES['Cython'] = 'cython'
    PYTEST_HEADER_MODULES['Scikit-image'] = 'skimage'
    PYTEST_HEADER_MODULES['asdf'] = 'asdf'
    PYTEST_HEADER_MODULES['pyarrow'] = 'pyarrow'
    TESTED_VERSIONS['Astropy'] = __version__


# This has to be in the root dir or it will not display in CI.
def pytest_report_header(config):
    # This gets added after the pytest-astropy-header output.
    return (f'ARCH_ON_CI: {os.environ.get("ARCH_ON_CI", "undefined")}\n'
            f'IS_CRON: {os.environ.get("IS_CRON", "undefined")}\n')


# Tell Hypothesis that we might be running slow tests, to print the seed blob
# so we can easily reproduce failures from CI, and derive a fuzzing profile
# to try many more inputs when we detect a scheduled build or when specifically
# requested using the HYPOTHESIS_PROFILE=fuzz environment variable or
# `pytest --hypothesis-profile=fuzz ...` argument.

hypothesis.settings.register_profile(
    'ci', deadline=None, print_blob=True, derandomize=True
)
hypothesis.settings.register_profile(
    'fuzzing', deadline=None, print_blob=True, max_examples=1000
)
default = 'fuzzing' if (os.environ.get('IS_CRON') == 'true' and os.environ.get('ARCH_ON_CI') not in ('aarch64', 'ppc64le')) else 'ci'  # noqa: E501
hypothesis.settings.load_profile(os.environ.get('HYPOTHESIS_PROFILE', default))

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
