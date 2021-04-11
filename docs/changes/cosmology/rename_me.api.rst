Directly instantiating the ``Cosmology`` base class now clones the current
default cosmology (controlled by ``default_cosmology``), updating parameters,
as specified. This allows for:

>>> from astropy.cosmology import Cosmology, default_cosmology
>>> print(Cosmology())  # the default cosmology
FlatLambdaCDM(name="Planck18", H0=67.7 km / (Mpc s), ...)

>>> with default_cosmology.set("WMAP5"):
>>>      print(Cosmology())  # gets the new default
FlatLambdaCDM(name="WMAP5", H0=70.2 km / (Mpc s), ...)

>>> Cosmology(name="modified", H0=70)  # perturb around the default
FlatLambdaCDM(name="modified", H0=70 km / (Mpc s), ...)