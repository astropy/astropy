# Licensed under a 3-clause BSD style license - see LICENSE.rst

from sys import version_info

if version_info[0]>2:
    from ..configobj_py3 import configobj,validate,__doc__
else:
    from ..configobj_py2 import configobj,validate,__doc__

del version_info #cleans up the namespace
