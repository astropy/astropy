============================
Writing Command-Line Scripts
============================

Command-line scripts in Astropy should follow a consistent scheme to promote
readability and compatibility.

Setuptools' `"entry points"`_ are used to automatically generate wrappers with
the correct extension. The scripts can live in their own module, or be part of
a larger module that implements a class or function for astropy library use.
They should have a ``main`` function to parse the arguments and pass those
arguments on to some library function so that the library function can be used
programatically when needed. The ``main`` function should accept an optional
single argument that holds the ``sys.argv`` list, except for the script name
(e.g., ``argv[1:]``). It must then be added to the list of entry points in the
``setup.py`` file (see the example below).

Command-line options can be parsed however desired, but the :mod:`argparse`
module is recommended when possible, due to its simpler and more flexible
interface relative to the older :mod:`optparse`.  :mod:`argparse` is only
available in python >=2.7 and >=3.2, however, so it should be imported as ``from
astropy.util.compat import argparse`` .

.. _"entry points": https://pythonhosted.org/setuptools/setuptools.html#automatic-script-creation

Example
-------

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

Then add the script to the ``setup.py`` ::

    entry_points['console_scripts'] = [
        'somescript = astropy.somepackage.somemod:main',
        ...
    ]
