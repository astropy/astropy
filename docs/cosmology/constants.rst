.. _astropy-cosmology-constants:

***************************
Cosmology-related Constants
***************************

.. currentmodule:: astropy.cosmology.constants

.. versionadded:: 6.0

The :mod:`astropy.cosmology.constants` contains constants used in the calculations in
the :mod:`astropy.cosmology` library. Currently the following constants are defined:


===========  ============================== ==================
Name         Description                    Units
===========  ============================== ==================
``c``        speed of light                 km s-1
``G``        gravitational constant         pc km2 s-2 Msun-1
===========  ============================== ==================


Cosmology API
-------------

Astropy cosmology is adopting the `Cosmology API`_, which includes adding this module
for constants. As more constants are added to the API,
:mod:`astropy.cosmology.constants` will be updated to include them.


Reference/API
=============

.. automodapi:: astropy.cosmology.constants
   :inherited-members:
