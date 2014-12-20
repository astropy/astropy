============================
Writing Command-Line Scripts
============================

Command-line scripts in Astropy should follow a consistent scheme to promote
readability and compatibility.

The actual script should be in the ``/scripts`` directory of the Astropy
source distribution, and should do nothing aside from importing a ``main``
function from astropy and execute it.  This was partly necessary because the
"2to3" utility that converted python 2.x code to 3.x does not convert scripts.
These scripts should be executable, include ``#!/usr/bin/env python`` at the
top, and should *not* end in ``.py``.

The ``main`` functions these scripts call should accept an optional single
argument that holds the ``sys.argv`` list, except for the script name
(e.g., ``argv[1:]``). This function can live in its own module, or be part of a
larger module that implements a class or function for astropy library use. The
``main`` function should do very little actual work - it should only parse the
arguments and pass those arguments on to some library function so that the
library function can be used programmatically when needed.
Command-line options can be parsed however desired, but the :mod:`argparse`
module is recommended when possible, due to its simpler and more flexible
interface relative to the older :mod:`optparse`. :mod:`argparse` is only
available in python >=2.7 and >=3.2, however, so it should be imported as
``from astropy.util.compat import argparse`` .


Example
-------

Contents of ``/scripts/cmdlinescript`` ::

    #!/usr/bin/env python
    # -*- coding: utf-8 -*-
    """An astropy command-line script"""

    import astropy.somepackage.somemod

    astropy.somepackage.somemod.main()

Contents of ``/astropy/somepackage/somemod.py`` ::

    def do_something(args, option=False):
        for a in args:
            if option:
                ...do something...
            else:
                ...do something else...

    def main(args=None):
        from astropy.utils.compat import argparse

        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument('-o', '--option', dest='op',action='store_true',
                            help='Some option that turns something on.')
        parser.add_argument('stuff', metavar='S', nargs='+',
                            help='Some input I should be able to get lots of.')

        res = parser.parse_args(args)

        do_something(res.stuff,res.op)

