:orphan:

.. _vision:

============================================
Vision for a Common Astronomy Python Package
============================================

The following document summarizes a vision for a common Astronomy Python
package, and how we can best all work together to achieve this. In the
following document, this common package will be referred to as the core
package. This vision is not set in stone, and we are committed to adapting it
to whatever process and guidelines work in practice.

The ultimate goal that we seek is a package that would contain much of the core
functionality and some common tools required across Astronomy, but not
*everything* Astronomers will ever need. The aim is primarily to avoid
duplication for common core tasks, and to provide a robust framework upon which
to build more complex tools.

Such a common package should not preclude any other Astronomy package from
existing, because there will always be more complex and/or specialized tools
required. These tools will be able to rely on a single core library for many
tasks, and thus reduce the number of dependencies, reduce duplication of
functionality, and increase consistency of their interfaces.

Procedure
---------

With the help of the community, the coordination committee will start by
identifying a few of key areas where initial development/consolidation will be
needed (such as FITS, WCS, coordinates, tables, photometry, spectra, etc.) and
will encourage teams to be formed to build standalone packages implementing
this functionality. These packages will be referred to as affiliated packages
(meaning that they are intended for future integration in the core package).

A set of requirements will be set out concerning the interfaces and
classes/methods that affiliated packages will need to make available in order
to ensure consistency between the different components. As the core package
grows, new potential areas/components for the core package will be identified.
Competition cannot be avoided, and will not be actively discouraged, but
whenever possible, developers should strive to work as a team to provide a
single and robust affiliated package, for the benefit of the community.

The affiliated packages will be developed outside the core package in
independent repositories, which will allow the teams the choice of tool and
organization. Once an affiliated package has implemented the desired
functionality, and satisfies quality criteria for coding style, documentation,
and testing, it will be considered for inclusion in the core package, and
further development will be done directly in the core package either via direct
access to the repository, or via patches/pull requests (exactly how this will
be done will be decided later).

To ensure uniformity across affiliated packages, and to facilitate integration
with the core package, developers who wish to submit their affiliated packages
for inclusion in the core will need to follow the layout of a ‘template’
package that will be provided before development starts.

Dependencies
------------

Affiliated packages should be able to be imported with only the following
dependencies:

* The Python Standard Library NumPy, SciPy, and Matplotlib Components already
  * in the core Astronomy package

Other packages may be used, but must be imported as needed rather than during
the initial import of the package.

If a dependency is needed, but is an affiliated package, the dependent package
will need to wait until the dependency is integrated into the core package
before being itself considered for inclusion. In the mean time, it can make use
of the other affiliated package in its current form, or other packages, so as
not to stall development. Thus, the first packages to be included in the core
will be those only requiring the standard library, NumPy, SciPy, and
Matplotlib.

If the required dependency will never be part of a main package, then by
default the dependency can be included but should be imported as needed
(meaning that it only prevents the importing of that component, not the entire
core package), unless a strong case is made and a general consensus is reached
by the community that this dependency is important enough to be required at a
higher level.

This system means that packages will be integrated into the core package in an
order depending on the dependency tree, and also ensures that the interfaces of
packages being integrated into the core package are consistent with those
already in the core package.

Initially, no dependency on GUI toolkits will be allowed in the core package.
If the community reaches agrees on a single toolkit that could be used, then
this toolkit will be allowed (but will only be imported as needed).

Keeping track of affiliated packages
------------------------------------

Affiliated packages will be listed in a central location (in addition to PyPI)
that will allow an easy installation of all the affiliated packages, for
example with a script that will seamlessly download and install all the
affiliated packages. The core package will also include mechanisms to
facilitate this installation process.

Existing Packages
-----------------

Developers who already have existing packages will be encouraged to continue
supporting them for the benefit of users until the core library is considered
stable, contains this functionality, and is released to the community.
Thereafter, developers should encourage users to transition to using the
functionality in the core package, and eventually phase out their own packages,
unless they provide added value over the core package.
