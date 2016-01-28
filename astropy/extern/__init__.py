# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This packages contains python packages that are bundled with Astropy but are
external to Astropy, and hence are developed in a separate source tree.  Note
that this package is distinct from the /cextern directory of the source code
distribution, as that directory only contains C extension code.

See the README.rst in this directory of the Astropy source repository for more
details.
"""

import sys
from distutils.version import StrictVersion as V


_VENDORED = [('six', '1.7.3')]
_SEARCH_PATH = ['astropy.extern.bundled.', '']


class VendorImporter(object):
    """
    A PEP 302 meta path importer for finding optionally-vendored
    or otherwise naturally-installed packages from __name__.
    """

    @staticmethod
    def _vendored_match(target):
        """
        Check if the given module/package name is in the _VENDORED
        module list, or is a subpackage of a module in the _VENDORED
        list.

        Returns the vendored module name and, if given, a required version
        for that module.
        """

        for v in _VENDORED:
            if isinstance(v, tuple):
                modname = v[0]
            else:
                modname = v

            if target == modname or target.startswith(modname + '.'):
                return v

    def find_module(self, fullname, path=None):
        root, base, target = fullname.partition(__name__ + '.')
        if root:
            return

        if not self._vendored_match(target):
            return

        return self

    def load_module(self, fullname):
        root, base, target = fullname.partition(__name__ + '.')
        vendored = self._vendored_match(target)
        version = None
        if isinstance(vendored, tuple):
            vendored, version = vendored

        for prefix in _SEARCH_PATH:
            try:
                __import__(prefix + target)
            except ImportError:
                continue

            mod = sys.modules[prefix + target]

            if target == vendored and version:
                # Only for top-level modules
                try:
                    if V(mod.__version__) >= V(version):
                        break
                except (AttributeError, ValueError):
                    # Attribute error if the module isn't what it should be and
                    # doesn't have a .__version__; ValueError if the version
                    # string exists but is somehow bogus/unparseable
                    continue
            else:
                break
        else:
            at_least = 'of minimum version {0} '.format(version)
            raise ImportError(
                "The '{target}' package {at_least}is required; "
                "normally this is bundled with this package so if you get "
                "this warning, consult the packager of your "
                "distribution.".format(**locals())
            )

        sys.modules[fullname] = mod
        return mod

    @classmethod
    def install(cls):
        if not any(isinstance(imp, cls) for imp in sys.meta_path):
            sys.meta_path.append(cls())


VendorImporter.install()
