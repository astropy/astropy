import sys

from ..hub_script import hub_script

from ..utils import ALLOW_INTERNET
ALLOW_INTERNET.set(False)


def test_hub_script():
    sys.argv.append('-m')  # run in multiple mode
    sys.argv.append('-w')  # disable web profile
    hub_script(timeout=3)
