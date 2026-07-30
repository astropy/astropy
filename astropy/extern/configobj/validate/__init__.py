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
from configobj.validate import *

