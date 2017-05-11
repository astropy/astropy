****************************
Getting started with subpackages
****************************

Because different subpackages have very different functionality, each subpackage has its own
getting started guide. These can be found by browsing the sections listed in the :ref:`user-docs`.

You can also look at docstrings for a
particular package or object, or access their documentation using the
`~astropy.utils.misc.find_api_page` function. For example, ::

    >>> from astropy import find_api_page
    >>> from astropy.units import Quantity
    >>> find_api_page(Quantity)  # doctest: +SKIP

will bring up the documentation for the `~astropy.units.Quantity` class
in your browser.
