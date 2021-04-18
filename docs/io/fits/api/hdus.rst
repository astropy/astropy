.. currentmodule:: astropy.io.fits

Header Data Unit
****************

Header Data Units are the fundamental container structure of the FITS format
consisting of a ``data`` member and its associated metadata in a ``header``.
They are defined in ``astropy.io.fits.hdu``.

The :class:`ImageHDU` and :class:`CompImageHDU` classes are discussed in the
section on :ref:`Images`.

The :class:`TableHDU` and :class:`BinTableHDU` classes are discussed in the
section on :ref:`Tables`.

:class:`PrimaryHDU`
===================
.. autoclass:: PrimaryHDU
   :members:
   :inherited-members:
   :show-inheritance:

:class:`GroupsHDU`
==================
.. autoclass:: GroupsHDU
   :members:
   :inherited-members:
   :show-inheritance:

:class:`GroupData`
==================
.. autoclass:: GroupData
   :members:
   :show-inheritance:

:class:`Group`
--------------
.. autoclass:: Group
   :members:
   :show-inheritance:

:class:`StreamingHDU`
=====================
.. autoclass:: StreamingHDU
   :members:
   :inherited-members:
   :show-inheritance:
