"""
This is a backwards compatibility-shim to support:

 ```
 import validate
 ```

 in a future release, we'd expect this to no longer work and
 instead using:

 ```
 import configobj.validate
 ```

 or:

 ```
 from configobj import validate
 ```
"""
# ASTROPY PATCH: relative import since this is vendored, not pip-installed
from ..configobj.validate import *

