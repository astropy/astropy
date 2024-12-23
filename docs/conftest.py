# Licensed under a 3-clause BSD style license - see LICENSE.rst

# This file needs to be included here to make sure commands such
# as ``pytest docs/...`` works, since this
# will ignore the conftest.py file at the root of the repository
# and the one in astropy/conftest.py

import os
from pathlib import Path

import pytest


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
