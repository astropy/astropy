# STDLIB
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))  # allows import of "mypackage"

# isort split
# THIRD PARTY
import mypackage

# LOCAL
from . import helper
