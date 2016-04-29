Developer guide to astropy.modeling
===================================

The purpose of this document is to describe the layout of the
astropy.modeling package as well as design choices, in order to make it
easier for people to contribute to it.

Code layout
-----------

The code in astropy.modeling is split into several categories:

Core code defining the base parameter and model classes:

* ``core.py``
* ``parameters.py``

Built-in models:

* ``functional_models.py``
* ``polynomial.py``
* ``powerlaws.py``
* ``projections.py``
* ``rotations.py``
* ``mappings.py``

Models module exposing all the model functionality above:

* ``models.py``

Code relating to fitting:

* ``fitting.py``
* ``optimizers.py``
* ``statistic.py``

General utilities:

* ``utils.py``
