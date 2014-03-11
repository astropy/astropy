import sys

if ((sys.version_info[0] == 2 and
     sys.version_info[1] >= 7) or
    (sys.version_info[0] == 3 and
     sys.version_info[1] >= 1)):
    from importlib import *
else:
    from ._importlib import *
