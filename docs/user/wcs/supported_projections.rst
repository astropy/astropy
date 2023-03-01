.. include:: references.txt

.. supported_projections:

Supported projections
---------------------

As `astropy.wcs` is based on `wcslib`_, it supports the standard
projections defined in the `FITS WCS standard`_.  These projection
codes are three letter strings specified in the second part of the ``CTYPEn`` keywords
(accessible through `Wcsprm.ctype <astropy.wcs.Wcsprm.ctype>`). For
example, a tangent projection with RA, DEC coordinates is defined by
``CTYPE1 = RA---TAN`` and ``CTYPE2 = DEC--TAN``. If a SIP distortion is present the
keywords become ``CTYPE1 = RA---TAN-SIP`` and ``CTYPE2 = DEC--TAN-SIP``.

The supported projection codes are:

- ``AZP``: zenithal/azimuthal perspective
- ``SZP``: slant zenithal perspective
- ``TAN``: gnomonic
- ``STG``: stereographic
- ``SIN``: orthographic/synthesis
- ``ARC``: zenithal/azimuthal equidistant
- ``ZPN``: zenithal/azimuthal polynomial
- ``ZEA``: zenithal/azimuthal equal area
- ``AIR``: Airy's projection
- ``CYP``: cylindrical perspective
- ``CEA``: cylindrical equal area
- ``CAR``: plate carrée
- ``MER``: Mercator's projection
- ``COP``: conic perspective
- ``COE``: conic equal area
- ``COD``: conic equidistant
- ``COO``: conic orthomorphic
- ``SFL``: Sanson-Flamsteed ("global sinusoid")
- ``PAR``: parabolic
- ``MOL``: Mollweide's projection
- ``AIT``: Hammer-Aitoff
- ``BON``: Bonne's projection
- ``PCO``: polyconic
- ``TSC``: tangential spherical cube
- ``CSC``: COBE quadrilateralized spherical cube
- ``QSC``: quadrilateralized spherical cube
- ``HPX``: HEALPix
- ``XPH``: HEALPix polar, aka "butterfly"

And, if built with wcslib 5.0 or later, the following polynomial
distortions are supported:

- ``TPV``: Polynomial distortion
- ``TUV``: Polynomial distortion

.. note::

    Though wcslib 5.4 and later handles ``SIP`` polynomial distortion,
    for backward compatibility, ``SIP`` is handled by astropy itself
    and methods exist to handle it specially.
