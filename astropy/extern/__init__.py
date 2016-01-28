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


_VENDORED = {'six': '1.7.3'}


class VendorImporter(object):
    """
    A PEP 302 meta path importer for finding optionally-vendored
    or otherwise naturally-installed packages from root_name.
    """

    def __init__(self, root_name, vendored_packages={}, vendor_pkg=None):
        self.root_name = root_name
        self.vendored_packages = vendored_packages
        self.vendor_pkg = vendor_pkg or root_name + 'extern.bundled'

    @property
    def search_path(self):
        """
        Search first the vendor package then as a natural package.
        """
        yield self.vendor_pkg + '.'
        yield ''

    def _vendored_match(self, target):
        """
        Check if the given module/package name is in the self.vendored_packages
        dict, or is a subpackage of a module self.vendored_packages.

        Returns the vendored module name and, if given, a required version for
        that module.
        """

        for v in self.vendored_packages:
            if target == v or target.startswith(v + '.'):
                return (v, self.vendored_packages[v])

    def find_module(self, fullname, path=None):
        root, base, target = fullname.partition(self.root_name + '.')
        if root:
            return

        if not self._vendored_match(target):
            return

        return self

    def load_module(self, fullname):
        root, base, target = fullname.partition(self.root_name + '.')
        vendored, version = self._vendored_match(target)

        for prefix in self.search_path:
            extant = prefix + target
            try:
                __import__(extant)
            except ImportError:
                continue

            mod = sys.modules[extant]

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
        # mysterious hack:
        # Remove the reference to the extant package/module
        # on later Python versions to cause relative imports
        # in the vendor package to resolve the same modules
        # as those going through this importer.
        if sys.version_info > (3, 3):
            del sys.modules[extant]
        return mod

    @classmethod
    def install(cls, *args, **kwargs):
        if not any(isinstance(imp, cls) for imp in sys.meta_path):
            sys.meta_path.append(cls(*args, **kwargs))


VendorImporter.install(__name__, _VENDORED)
