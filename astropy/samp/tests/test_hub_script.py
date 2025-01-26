import logging
import sys

import pytest

from astropy.samp import conf
from astropy.samp.hub_script import hub_script


def setup_module(module):
    conf.use_internet = False


def setup_function(function):
    function.sys_argv_orig = sys.argv
    sys.argv = ["samp_hub"]


def teardown_function(function):
    sys.argv = function.sys_argv_orig


@pytest.mark.slow
def test_hub_script(monkeypatch):
    mock_logger = logging.getLogger(__name__)
    monkeypatch.setattr("astropy.samp.hub_script.log", mock_logger)
    sys.argv.append("-m")  # run in multiple mode
    sys.argv.append("-w")  # disable web profile
    sys.argv.extend(["--log-level", "CRITICAL"])  # set the logging level

    hub_script(timeout=3)

    assert mock_logger.level == logging.getLevelNamesMapping()["CRITICAL"]
