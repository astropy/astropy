# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file needs to be included here to make sure commands such
# as ``pytest docs/...`` works, since this
# will ignore the conftest.py file at the root of the repository
# and the one in astropy/conftest.py

import os
from pathlib import Path

import pytest

# Make sure we use temporary directories for the config and cache
# so that the tests are insensitive to local configuration. Note that this
# is also set in the test runner, but we need to also set it here for
# things to work properly in parallel mode
# note: session-level + autouse doesn't require a cleanup phase


@pytest.fixture(scope="session", autouse=True)
def _session_level_cache_dir(tmp_path_factory):
    os.environ["ASTROPY_CACHE_DIR"] = str(tmp_path_factory.mktemp("astropy_cache_"))
    os.environ["XDG_CACHE_HOME"] = str(tmp_path_factory.mktemp("xdg_cache_"))


@pytest.fixture(scope="session", autouse=True)
def _session_level_config_dir(tmp_path_factory):
    os.environ["ASTROPY_CONFIG_DIR"] = str(tmp_path_factory.mktemp("astropy_config_"))
    os.environ["XDG_CONFIG_HOME"] = str(tmp_path_factory.mktemp("xdg_config_"))


@pytest.fixture(autouse=True)
def _docdir(request):
    """Run doctests in isolated tmp_path so outputs do not end up in repo."""
    # Trigger ONLY for doctestplus
    doctest_plugin = request.config.pluginmanager.getplugin("doctestplus")
    if isinstance(request.node.parent, doctest_plugin._doctest_textfile_item_cls):
        # Don't apply this fixture to io.rst.  It reads files and doesn't write.
        # Implementation from https://github.com/pytest-dev/pytest/discussions/10437
        if "io.rst" not in request.node.name:
            old_cwd = Path.cwd()
            tmp_path = request.getfixturevalue("tmp_path")
            os.chdir(tmp_path)
            yield
            os.chdir(old_cwd)
        else:
            yield
    else:
        yield
