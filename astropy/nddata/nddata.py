# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This module implements the base NDData class.


from copy import deepcopy

import numpy as np

from astropy import log
from astropy.units import Quantity, Unit
from astropy.utils.compat.optional_deps import HAS_DASK
from astropy.utils.masked import Masked, MaskedNDArray
from astropy.utils.metadata import MetaData
from astropy.wcs.wcsapi import (
    BaseHighLevelWCS,
    BaseLowLevelWCS,
    HighLevelWCSWrapper,
    SlicedLowLevelWCS,  # noqa: F401
)

from .nddata_base import NDDataBase
from .nduncertainty import NDUncertainty, UnknownUncertainty

__all__ = ["NDData"]

_meta_doc = """`dict`-like : Additional meta information about the dataset."""


class NDData(NDDataBase):
    """
    A container for `numpy.ndarray`-based datasets, using the
    `~astropy.nddata.NDDataBase` interface.

    The key distinction from raw `numpy.ndarray` is the presence of
    additional metadata such as uncertainty, mask, unit, a coordinate system
    and/or a dictionary containing further meta information. This class *only*
    provides a container for *storing* such datasets. For further functionality
    take a look at the ``See also`` section.

    See also: https://docs.astropy.org/en/stable/nddata/

    Parameters
    ----------
    data : `numpy.ndarray`-like or `NDData`-like
        The dataset.

    uncertainty : any type, optional
        Uncertainty in the dataset.
        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, for example ``"std"`` for standard deviation or
        ``"var"`` for variance. A metaclass defining such an interface is
        `NDUncertainty` - but isn't mandatory. If the uncertainty has no such
        attribute the uncertainty is stored as `UnknownUncertainty`.
        Defaults to ``None``.

    mask : any type, optional
        Mask for the dataset. Masks should follow the ``numpy`` convention that
        **valid** data points are marked by ``False`` and **invalid** ones with
        ``True``.
        Defaults to ``None``.

    wcs : any type, optional
        World coordinate system (WCS) for the dataset.
        Default is ``None``.

    meta : `dict`-like object, optional
        Additional meta information about the dataset. If no meta is provided
        an empty dict is created.
        Default is ``None``.

    unit : unit-like, optional
        Unit for the dataset. Strings that can be converted to a
        `~astropy.units.Unit` are allowed.
        Default is ``None``.

    copy : `bool`, optional
        Indicates whether to save the arguments as copy. ``True`` copies
        every attribute before saving it while ``False`` tries to save every
        parameter as reference.
        Note however that it is not always possible to save the input as
        reference.
        Default is ``False``.

        .. versionadded:: 1.2

    psf : `numpy.ndarray` or None, optional
        Image representation of the PSF. In order for convolution to be flux-
        preserving, this should generally be normalized to sum to unity.

    Raises
    ------
    TypeError
        In case ``data`` or ``meta`` don't meet the restrictions.

    Notes
    -----
    Each attribute can be accessed through the homonymous instance attribute:
    ``data`` in a `NDData` object can be accessed through the `data`
    attribute::

        >>> from astropy.nddata import NDData
        >>> nd = NDData([1,2,3])
        >>> nd.data
        array([1, 2, 3])

    Given a conflicting implicit and an explicit parameter during
    initialization, for example the ``data`` is a `~astropy.units.Quantity` and
    the unit parameter is not ``None``, then the implicit parameter is replaced
    (without conversion) by the explicit one and a warning is issued::

        >>> import numpy as np
        >>> import astropy.units as u
        >>> q = np.array([1,2,3,4]) * u.m
        >>> nd2 = NDData(q, unit=u.cm)
        INFO: overwriting Quantity's current unit with specified unit. [astropy.nddata.nddata]
        >>> nd2.data  # doctest: +FLOAT_CMP
        array([100., 200., 300., 400.])
        >>> nd2.unit
        Unit("cm")

    See Also
    --------
    NDDataRef
    NDDataArray
    """

    # Instead of a custom property use the MetaData descriptor also used for
    # Tables. It will check if the meta is dict-like or raise an exception.
    meta = MetaData(doc=_meta_doc, copy=False)

    def __init__(
        self,
        data,
        uncertainty=None,
        mask=None,
        wcs=None,
        meta=None,
        unit=None,
        copy=False,
        psf=None,
    ):
        # Rather pointless since the NDDataBase does not implement any setting
        # but before the NDDataBase did call the uncertainty
        # setter. But if anyone wants to alter this behavior again the call
        # to the superclass NDDataBase should be in here.
        super().__init__()

        # Check if data is any type from which to collect some implicitly
        # passed parameters.
        if isinstance(data, NDData):  # don't use self.__class__ (issue #4137)
            # Of course we need to check the data because subclasses with other
            # init-logic might be passed in here. We could skip these
            # tests if we compared for self.__class__ but that has other
            # drawbacks.

            # Comparing if there is an explicit and an implicit unit parameter.
            # If that is the case use the explicit one and issue a warning
            # that there might be a conflict. In case there is no explicit
            # unit just overwrite the unit parameter with the NDData.unit
            # and proceed as if that one was given as parameter. Same for the
            # other parameters.
            if unit is None and data.unit is not None:
                unit = data.unit
            elif unit is not None and data.unit is not None:
                log.info("overwriting NDData's current unit with specified unit.")

            if uncertainty is not None and data.uncertainty is not None:
                log.info(
                    "overwriting NDData's current "
                    "uncertainty with specified uncertainty."
                )
            elif data.uncertainty is not None:
                uncertainty = data.uncertainty

            if mask is not None and data.mask is not None:
                log.info("overwriting NDData's current mask with specified mask.")
            elif data.mask is not None:
                mask = data.mask

            if wcs is not None and data.wcs is not None:
                log.info("overwriting NDData's current wcs with specified wcs.")
            elif data.wcs is not None:
                wcs = data.wcs

            if psf is not None and data.psf is not None:
                log.info("Overwriting NDData's current psf with specified psf.")
            elif data.psf is not None:
                psf = data.psf

            if meta is not None and data.meta is not None:
                log.info("overwriting NDData's current meta with specified meta.")
            elif data.meta is not None:
                meta = data.meta

            # get the data attribute as it is, and continue to process it:
            data = data.data

        # if the data is wrapped by astropy.utils.masked.Masked:
        if isinstance(data, Masked):
            # first get the mask if one is available:
            if hasattr(data, "mask"):
                if mask is not None:
                    log.info(
                        "overwriting Masked Quantity's current mask with specified mask."
                    )
                else:
                    mask = data.mask

            if isinstance(data, MaskedNDArray):
                if unit is not None and hasattr(data, "unit") and data.unit != unit:
                    log.info(
                        "overwriting MaskedNDArray's current unit with specified unit."
                    )
                    data = data.to(unit).value
                elif unit is None and hasattr(data, "unit"):
                    unit = data.unit
                    data = data.value

                # now get the unmasked ndarray:
                data = np.asarray(data)

            if isinstance(data, Quantity):
                # this is a Quantity:
                if unit is not None and data.unit != unit:
                    log.info("overwriting Quantity's current unit with specified unit.")
                    data = data.to(unit)
                elif unit is None and data.unit is not None:
                    unit = data.unit
                data = data.value

        if isinstance(data, np.ma.masked_array):
            if mask is not None:
                log.info(
                    "overwriting masked ndarray's current mask with specified mask."
                )
            else:
                mask = data.mask
            data = data.data

        if isinstance(data, Quantity):
            # this is a Quantity:
            if unit is not None and data.unit != unit:
                log.info("overwriting Quantity's current unit with specified unit.")
                data = data.to(unit)
            elif unit is None and data.unit is not None:
                unit = data.unit
            data = data.value

        if isinstance(data, np.ndarray):
            # check for mask from np.ma.masked_ndarray
            if hasattr(data, "mask"):
                if mask is not None:
                    log.info(
                        "overwriting masked ndarray's current mask with specified mask."
                    )
                else:
                    mask = data.mask

        # Quick check on the parameters if they match the requirements.
        if (
            not hasattr(data, "shape")
            or not hasattr(data, "__getitem__")
            or not hasattr(data, "__array__")
        ):
            # Data doesn't look like a numpy array, try converting it to
            # one.
            data = np.asanyarray(data)
        # Another quick check to see if what we got looks like an array
        # rather than an object (since numpy will convert a
        # non-numerical/non-string inputs to an array of objects).
        if data.dtype == "O":
            raise TypeError("could not convert data to numpy array.")

        if unit is not None:
            unit = Unit(unit)

        if copy:
            # Data might have been copied before but no way of validating
            # without another variable.
            data = deepcopy(data)
            mask = deepcopy(mask)
            wcs = deepcopy(wcs)
            psf = deepcopy(psf)
            meta = deepcopy(meta)
            uncertainty = deepcopy(uncertainty)
            # Actually - copying the unit is unnecessary but better safe
            # than sorry :-)
            unit = deepcopy(unit)

        # Store the attributes
        self._data = data
        self.mask = mask
        self._wcs = None
        if wcs is not None:
            # Validate the wcs
            self.wcs = wcs
        self.meta = meta  # TODO: Make this call the setter sometime
        self._unit = unit
        # Call the setter for uncertainty to further check the uncertainty
        self.uncertainty = uncertainty
        self.psf = psf

    def __str__(self):
        data = str(self.data)
        unit = f" {self.unit}" if self.unit is not None else ""

        return data + unit

    def __repr__(self):
        prefix = self.__class__.__name__ + "("
        # to support reprs for other non-ndarray `data` attributes,
        # add more cases here:
        if HAS_DASK:
            import dask.array as da

            is_dask = isinstance(self.data, da.Array)
        else:
            is_dask = False

        if (
            isinstance(self.data, (int, float, np.ndarray))
            or np.issubdtype(float, self.data)
            or np.issubdtype(int, self.data)
        ) and not is_dask:
            # if data is an ndarray, get build a repr via Masked:
            ma = Masked(self.data, mask=self.mask)
            data_repr = repr(ma)
            # strip the class name generated by Masked, and the final ")",
            # and re-pad multi-line array outputs
            after_first_paren = data_repr.index("(") + 1
            old_prefix_spaces = " " * after_first_paren
            new_prefix_spaces = " " * len(prefix)
            data_repr = data_repr[after_first_paren:-1].replace(
                old_prefix_spaces, new_prefix_spaces
            )
            unit = f", unit='{self.unit}'" if self.unit is not None else ""
            return f"{prefix}{data_repr}{unit})"
        else:
            # relying on the builtin repr for objects other than ndarrays
            # allows for compatibility with e.g. dask:
            contents = []
            for attr in ("data", "mask", "uncertainty", "unit"):
                attr_data = getattr(self, attr)
                if attr_data is not None:
                    attr_prefix = f"\n  {attr}="
                    attr_repr = repr(attr_data)
                    attr_repr = attr_repr.replace(
                        "\n", f"\n{' ' * (len(attr_prefix) - 1)}"
                    )
                    contents.append(attr_prefix + attr_repr)
            return prefix + ",".join(contents) + "\n)"

    @property
    def data(self):
        """
        `~numpy.ndarray`-like : The stored dataset.
        """
        return self._data

    @property
    def mask(self):
        """
        any type : Mask for the dataset, if any.

        Masks should follow the ``numpy`` convention that valid data points are
        marked by ``False`` and invalid ones with ``True``.
        """
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def unit(self):
        """
        `~astropy.units.Unit` : Unit for the dataset, if any.
        """
        return self._unit

    @property
    def wcs(self):
        """
        any type : A world coordinate system (WCS) for the dataset, if any.
        """
        return self._wcs

    @wcs.setter
    def wcs(self, wcs):
        if self._wcs is not None and wcs is not None:
            raise ValueError(
                "You can only set the wcs attribute with a WCS if no WCS is present."
            )

        if wcs is None or isinstance(wcs, BaseHighLevelWCS):
            self._wcs = wcs
        elif isinstance(wcs, BaseLowLevelWCS):
            self._wcs = HighLevelWCSWrapper(wcs)
        else:
            raise TypeError(
                "The wcs argument must implement either the high or low level WCS API."
            )

    @property
    def psf(self):
        return self._psf

    @psf.setter
    def psf(self, value):
        self._psf = value

    @property
    def uncertainty(self):
        """
        any type : Uncertainty in the dataset, if any.

        Should have an attribute ``uncertainty_type`` that defines what kind of
        uncertainty is stored, such as ``'std'`` for standard deviation or
        ``'var'`` for variance. A metaclass defining such an interface is
        `~astropy.nddata.NDUncertainty` but isn't mandatory.
        """
        return self._uncertainty

    @uncertainty.setter
    def uncertainty(self, value):
        if value is not None:
            # There is one requirements on the uncertainty: That
            # it has an attribute 'uncertainty_type'.
            # If it does not match this requirement convert it to an unknown
            # uncertainty.
            if not hasattr(value, "uncertainty_type"):
                log.info("uncertainty should have attribute uncertainty_type.")
                value = UnknownUncertainty(value, copy=False)

            # If it is a subclass of NDUncertainty we must set the
            # parent_nddata attribute. (#4152)
            if isinstance(value, NDUncertainty):
                # In case the uncertainty already has a parent create a new
                # instance because we need to assume that we don't want to
                # steal the uncertainty from another NDData object
                if value._parent_nddata is not None:
                    value = value.__class__(value, copy=False)
                # Then link it to this NDData instance (internally this needs
                # to be saved as weakref but that's done by NDUncertainty
                # setter).
                value.parent_nddata = self
        self._uncertainty = value
