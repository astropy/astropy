# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from abc import ABCMeta, abstractproperty, abstractmethod
from copy import deepcopy

from ..utils.compat import ignored
from ..units import Unit, Quantity
from ..extern import six
from .. import log

__all__ = ['MissingDataAssociationException',
           'IncompatibleUncertaintiesException', 'NDUncertainty',
           'StdDevUncertainty']


class IncompatibleUncertaintiesException(Exception):
    """
    This exception should be used to indicate cases in which uncertainties
    with two different classes can not be propagated.
    """


class MissingDataAssociationException(Exception):
    """
    This exception should be used to indicate that an uncertainty instance has
    not been associated with a parent `~astropy.nddata.NDData` object.
    """


@six.add_metaclass(ABCMeta)
class NDUncertainty(object):
    """
    This is the base class for uncertainty classes used with NDData.

    Parameters
    ----------
    array: any type, optional
        The array or value (the parameter name is due to historical reasons) of
        the uncertainty. There is nothing enforced but `~numpy.ndarray`,
        `~astropy.units.Quantity` or `NDUncertainty` subclasses are
        recommended.
    unit: `~astropy.units.Unit` or `str`, optional
        The unit of the uncertainty `array`. If input is a `str` this will be
        converted to an `~astropy.units.Unit`.
    copy: `bool`, optional
        Should the uncertainty be saved as reference (``False``) or as a copy
        (``True``). Defaults to True.

    Notes
    -----
    1. NDUncertainty is an abstract class and should *never* be instantiated
       directly.

    2. NDUncertainty provides a ``__getitem__`` method so slicing is in theory
       supported but that relies on the `array` being something that actually
       can be sliced.

    3. Subclasses need to define:

       - property ``uncertainty_type`` which should be a small string defining
         the kind of uncertainty. Recommened is ``std`` for standard deviation
         ``var`` for variance (following the ``numpy`` conventions).
       - tbc...

    4. NDUncertainty and it's subclasses can only be used for uncertainty
       propagation if their ``parent_nddata`` is a reference to the `NDData`
       object. This is needed since many kinds of propagation need the actual
       data besides the uncertainty.
    """

    # Indicates whether the class supports the propagation of correlated
    # uncertainties
    supports_correlated = False

    def __init__(self, array=None, unit=None, copy=True):
        # It is possible to create an uncertainty from another Uncertainty or
        # a quantity
        if isinstance(array, StdDevUncertainty):
            if unit is not None and array.unit is not None:
                log.info("Overwriting Uncertainty's current "
                         "unit with specified unit")
            elif array.unit is not None:
                unit = array.unit
            array = array.array

        elif isinstance(array, Quantity):
            if unit is not None and array.unit is not None:
                log.info("Overwriting Quantity's current "
                         "unit with specified unit")
            elif array.unit is not None:
                unit = array.unit
            array = array.value

        if copy:
            array = deepcopy(array)
            unit = deepcopy(unit)

        self.array = array
        self.unit = unit

    def __getitem__(self, item):
        return self.__class__(self.array[item], unit=self.unit, copy=False)

    @abstractproperty
    def uncertainty_type(self):
        """
        `str`: Short description which kind of uncertainty is saved

        Defined as abstract property so subclasses *have* to override this
        and return a string.
        """
        return None

    @property
    def parent_nddata(self):
        """
        `NDData` reference: The `NDData` whose uncertainty this is.

        In case the reference is not set uncertainty propagation will not be
        possible since almost all kinds of propagation need the uncertain
        data besides the uncertainty.
        """
        message = "Uncertainty is not associated with an NDData object"
        try:
            if self._parent_nddata is None:
                raise MissingDataAssociationException(message)
            else:
                return self._parent_nddata
        except AttributeError:
            raise MissingDataAssociationException(message)

    @parent_nddata.setter
    def parent_nddata(self, value):
        self._parent_nddata = value

    @property
    def array(self):
        """
        any type: the uncertainty

        `numpy.ndarray`s or scalars are preferred since any kind of
        uncertainty propagation with other types is not defined in this meta
        class.
        """
        return self._array

    @array.setter
    def array(self, value):
        self._array = value

    @property
    def unit(self):
        """
        `~astropy.units.Unit`: The unit of the uncertainty

        Even though it is not enforced the unit should be convertable to the
        ``parent_nddata``s unit. Otherwise uncertainty propagation might give
        some wrong results. If the unit is not set the unit of the parent will
        be assumed.
        """
        if self._unit is None:
            if self.parent_nddata.unit is None:
                return None
            else:
                return self.parent_nddata.unit
        return self._unit

    @unit.setter
    def unit(self, value):
        if value is None:
            self._unit = None
        else:
            self._unit = Unit(value)

    @abstractmethod
    def propagate_add(self, other_nddata, result_data):
        return None

    @abstractmethod
    def propagate_subtract(self, other_nddata, result_data):
        return None

    @abstractmethod
    def propagate_multiply(self, other_nddata, result_data):
        return None

    @abstractmethod
    def propagate_divide(self, other_nddata, result_data):
        return None


class StdDevUncertainty(NDUncertainty):
    """
    A class for standard deviation uncertainties
    """

    support_correlated = False

    def __init__(self, array=None, unit=None, copy=True):
        self._unit = None
        if array is None:
            self.array = None
        elif isinstance(array, StdDevUncertainty):
            self.array = np.array(array.array, copy=copy, subok=True)
        elif isinstance(array, Quantity):
            self.array = np.array(array.value, copy=copy, subok=True)
            self._unit = array.unit
        else:
            self.array = np.array(array, copy=copy, subok=True)

    @property
    def uncertainty_type(self):
        return 'std'

    def propagate_add(self, other_nddata, result_data):
        """
        Propagate uncertainties for addition.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the
            uncertainties.
        """

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set "
                             "in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = np.sqrt(self.array**2 +
                                           other_nddata.uncertainty.array**2)

        return result_uncertainty

    def propagate_subtract(self, other_nddata, result_data):
        """
        Propagate uncertainties for subtraction.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the
            uncertainties.
        """

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set "
                             "in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = np.sqrt(self.array**2 +
                                           other_nddata.uncertainty.array**2)

        return result_uncertainty

    def propagate_multiply(self, other_nddata, result_data):
        """
        Propagate uncertainties for multiplication.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the
            uncertainties.
        """

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set in "
                             "other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = \
            (np.sqrt((self.array/self.parent_nddata.data)**2
             + (other_nddata.uncertainty.array/other_nddata.data)**2) *
             result_data)

        return result_uncertainty

    def propagate_divide(self, other_nddata, result_data):
        """
        Propagate uncertainties for division.

        Parameters
        ----------
        other_nddata : NDData instance
            The data for the second other_nddata in a + b
        result_data : `~numpy.ndarray` instance
            The data array that is the result of the addition

        Returns
        -------
        result_uncertainty : NDUncertainty instance
            The resulting uncertainty

        Raises
        ------
        IncompatibleUncertaintiesException
            Raised if the method does not know how to propagate the
            uncertainties.
        """

        if not isinstance(other_nddata.uncertainty, StdDevUncertainty):
            raise IncompatibleUncertaintiesException

        if self.array is None:
            raise ValueError("standard deviation values are not set")

        if other_nddata.uncertainty.array is None:
            raise ValueError("standard deviation values are not set "
                             "in other_nddata")

        result_uncertainty = StdDevUncertainty()
        result_uncertainty.array = \
            (np.sqrt((self.array/self.parent_nddata.data)**2
             + (other_nddata.uncertainty.array/other_nddata.data)**2) *
             result_data)

        return result_uncertainty
