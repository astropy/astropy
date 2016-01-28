"""Support WeakSet on Python 2.6"""

try:
    from weakref import WeakSet
except ImportError:
    from ._weakset_py2 import WeakSet
