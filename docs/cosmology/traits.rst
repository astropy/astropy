.. _astropy-cosmology-traits:

****************
Cosmology Traits
****************

.. currentmodule:: astropy.cosmology.traits

The ``traits`` module hosts various parts of cosmologies, such as the
:class:`~astropy.cosmology.traits.ScaleFactor`. These traits can be used to more easily
construct custom cosmologies by combining different components.

As a simple example, the :class:`~astropy.cosmology.traits.ScaleFactor` trait provides
the ``scale_factor0`` property and ``scale_factor0(z)`` method for computing the
cosmological scale factor, which is defined as :math:`a = a_0 / (1 + z)`. By using this
trait, you can add scale factor functionality to your custom cosmology class without
having to implement it from scratch.

Here is an example of how to use the :class:`~astropy.cosmology.traits.ScaleFactor`
trait in a custom cosmology class:

.. code-block:: python

    from astropy.cosmology.traits import ScaleFactor
    from astropy.cosmology import Cosmology

    class CustomCosmology(Cosmology, ScaleFactor):
        def __init__(self, H0, Om0, Ode0):
            self.H0 = H0
            self.Om0 = Om0
            self.Ode0 = Ode0
            super().__init__()

        # Additional custom methods and properties can be added here

    cosmo = CustomCosmology(H0=70, Om0=0.3, Ode0=0.7)
    print(cosmo.scale_factor(0))  # Output: 1.0

By combining different traits, you can create fully-featured cosmology classes with
minimal effort.


Reference/API
*************

.. py:currentmodule:: astropy.cosmology.traits

.. automodapi:: astropy.cosmology.traits
   :inherited-members:
