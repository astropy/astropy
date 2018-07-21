# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

# Dependencies
import numpy as np
import warnings

# Project
from .. import units as u
from ..utils.compat.numpy import broadcast_to as np_broadcast_to
from ..utils.exceptions import AstropyDeprecationWarning
from ..utils import OrderedDescriptor, ShapedLikeNDArray

__all__ = ['Attribute', 'TimeAttribute', 'QuantityAttribute',
           'EarthLocationAttribute', 'CoordinateAttribute',
           'CartesianRepresentationAttribute',
           'DifferentialAttribute']


class Attribute(OrderedDescriptor):
    """A non-mutable data descriptor to hold a frame attribute.

    This class must be used to define frame attributes (e.g. ``equinox`` or
    ``obstime``) that are included in a frame class definition.

    Examples
    --------
    The `~astropy.coordinates.FK4` class uses the following class attributes::

      class FK4(BaseCoordinateFrame):
          equinox = TimeAttribute(default=_EQUINOX_B1950)
          obstime = TimeAttribute(default=None,
                                  secondary_attribute='equinox')

    This means that ``equinox`` and ``obstime`` are available to be set as
    keyword arguments when creating an ``FK4`` class instance and are then
    accessible as instance attributes.  The instance value for the attribute
    must be stored in ``'_' + <attribute_name>`` by the frame ``__init__``
    method.

    Note in this example that ``equinox`` and ``obstime`` are time attributes
    and use the ``TimeAttributeFrame`` class.  This subclass overrides the
    ``convert_input`` method to validate and convert inputs into a ``Time``
    object.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    _class_attribute_ = 'frame_attributes'
    _name_attribute_ = 'name'
    name = '<unbound>'

    def __init__(self, default=None, secondary_attribute=''):
        self.default = default
        self.secondary_attribute = secondary_attribute
        super(Attribute, self).__init__()

    def convert_input(self, value):
        """
        Validate the input ``value`` and convert to expected attribute class.

        The base method here does nothing, but subclasses can implement this
        as needed.  The method should catch any internal exceptions and raise
        ValueError with an informative message.

        The method returns the validated input along with a boolean that
        indicates whether the input value was actually converted.  If the input
        value was already the correct type then the ``converted`` return value
        should be ``False``.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        output_value
            The ``value`` converted to the correct type (or just ``value`` if
            ``converted`` is False)
        converted : bool
            True if the conversion was actually performed, False otherwise.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        return value, False

    def __get__(self, instance, frame_cls=None):
        if instance is None:
            out = self.default
        else:
            out = getattr(instance, '_' + self.name, self.default)
            if out is None:
                out = getattr(instance, self.secondary_attribute, self.default)

        out, converted = self.convert_input(out)
        if instance is not None:
            instance_shape = getattr(instance, 'shape', None)
            if instance_shape is not None and (getattr(out, 'size', 1) > 1 and
                                               out.shape != instance_shape):
                # If the shapes do not match, try broadcasting.
                try:
                    if isinstance(out, ShapedLikeNDArray):
                        out = out._apply(np_broadcast_to, shape=instance_shape,
                                         subok=True)
                    else:
                        out = np_broadcast_to(out, instance_shape, subok=True)
                except ValueError:
                    # raise more informative exception.
                    raise ValueError(
                        "attribute {0} should be scalar or have shape {1}, "
                        "but is has shape {2} and could not be broadcast."
                        .format(self.name, instance_shape, out.shape))

                converted = True

            if converted:
                setattr(instance, '_' + self.name, out)

        return out

    def __set__(self, instance, val):
        raise AttributeError('Cannot set frame attribute')


class TimeAttribute(Attribute):
    """
    Frame attribute descriptor for quantities that are Time objects.
    See the `~astropy.coordinates.Attribute` API doc for further
    information.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def convert_input(self, value):
        """
        Convert input value to a Time object and validate by running through
        the Time constructor.  Also check that the input was a scalar.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """

        from ..time import Time

        if value is None:
            return None, False

        if isinstance(value, Time):
            out = value
            converted = False
        else:
            try:
                out = Time(value)
            except Exception as err:
                raise ValueError(
                    'Invalid time input {0}={1!r}\n{2}'.format(self.name,
                                                               value, err))
            converted = True

        return out, converted


class CartesianRepresentationAttribute(Attribute):
    """
    A frame attribute that is a CartesianRepresentation with specified units.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    unit : unit object or None
        Name of a unit that the input will be converted into. If None, no
        unit-checking or conversion is performed
    """

    def __init__(self, default=None, secondary_attribute='', unit=None):
        super(CartesianRepresentationAttribute, self).__init__(
            default, secondary_attribute)
        self.unit = unit

    def convert_input(self, value):
        """
        Checks that the input is a CartesianRepresentation with the correct
        unit, or the special value ``[0, 0, 0]``.

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """

        if (isinstance(value, list) and len(value) == 3 and
                all(v == 0 for v in value) and self.unit is not None):
            return CartesianRepresentation(np.zeros(3) * self.unit), True
        else:
            # is it a CartesianRepresentation with correct unit?
            if hasattr(value, 'xyz') and value.xyz.unit == self.unit:
                return value, False

            converted = True
            # if it's a CartesianRepresentation, get the xyz Quantity
            value = getattr(value, 'xyz', value)
            if not hasattr(value, 'unit'):
                raise TypeError('tried to set a {0} with something that does '
                                'not have a unit.'
                                .format(self.__class__.__name__))

            value = value.to(self.unit)

            # now try and make a CartesianRepresentation.
            cartrep = CartesianRepresentation(value, copy=False)
            return cartrep, converted


class QuantityAttribute(Attribute):
    """
    A frame attribute that is a quantity with specified units and shape
    (optionally).

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    unit : unit object or None
        Name of a unit that the input will be converted into. If None, no
        unit-checking or conversion is performed
    shape : tuple or None
        If given, specifies the shape the attribute must be
    """

    def __init__(self, default=None, secondary_attribute='', unit=None, shape=None):
        super(QuantityAttribute, self).__init__(default, secondary_attribute)
        self.unit = unit
        self.shape = shape

    def convert_input(self, value):
        """
        Checks that the input is a Quantity with the necessary units (or the
        special value ``0``).

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        if np.all(value == 0) and self.unit is not None:
            return u.Quantity(np.zeros(self.shape), self.unit), True
        else:
            if not hasattr(value, 'unit') and self.unit != u.dimensionless_unscaled:
                raise TypeError('Tried to set a QuantityAttribute with '
                                'something that does not have a unit.')
            oldvalue = value
            value = u.Quantity(oldvalue, self.unit, copy=False)
            if self.shape is not None and value.shape != self.shape:
                raise ValueError('The provided value has shape "{0}", but '
                                 'should have shape "{1}"'.format(value.shape,
                                                                  self.shape))
            converted = oldvalue is not value
            return value, converted

class EarthLocationAttribute(Attribute):
    """
    A frame attribute that can act as a `~astropy.coordinates.EarthLocation`.
    It can be created as anything that can be transformed to the
    `~astropy.coordinates.ITRS` frame, but always presents as an `EarthLocation`
    when accessed after creation.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def convert_input(self, value):
        """
        Checks that the input is a Quantity with the necessary units (or the
        special value ``0``).

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """

        if value is None:
            return None, False
        elif isinstance(value, EarthLocation):
            return value, False
        else:
            # we have to do the import here because of some tricky circular deps
            from .builtin_frames import ITRS

            if not hasattr(value, 'transform_to'):
                raise ValueError('"{0}" was passed into an '
                                 'EarthLocationAttribute, but it does not have '
                                 '"transform_to" method'.format(value))
            itrsobj = value.transform_to(ITRS)
            return itrsobj.earth_location, True


class CoordinateAttribute(Attribute):
    """
    A frame attribute which is a coordinate object. It can be given as a
    low-level frame class *or* a `~astropy.coordinates.SkyCoord`, but will
    always be converted to the low-level frame class when accessed.

    Parameters
    ----------
    frame : a coordinate frame class
        The type of frame this attribute can be
    default : object
        Default value for the attribute if not provided
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def __init__(self, frame, default=None, secondary_attribute=''):
        self._frame = frame
        super(CoordinateAttribute, self).__init__(default, secondary_attribute)

    def convert_input(self, value):
        """
        Checks that the input is a SkyCoord with the necessary units (or the
        special value ``None``).

        Parameters
        ----------
        value : object
            Input value to be converted.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """
        if value is None:
            return None, False
        elif isinstance(value, self._frame):
            return value, False
        else:
            if not hasattr(value, 'transform_to'):
                raise ValueError('"{0}" was passed into a '
                                 'CoordinateAttribute, but it does not have '
                                 '"transform_to" method'.format(value))
            transformedobj = value.transform_to(self._frame)
            if hasattr(transformedobj, 'frame'):
                transformedobj = transformedobj.frame
            return transformedobj, True


class DifferentialAttribute(Attribute):
    """A frame attribute which is a differential instance.

    The optional ``allowed_classes`` argument allows specifying a restricted
    set of valid differential classes to check the input against. Otherwise,
    any `~astropy.coordinates.BaseDifferential` subclass instance is valid.

    Parameters
    ----------
    default : object
        Default value for the attribute if not provided
    allowed_classes : tuple, optional
        A list of allowed differential classes for this attribute to have.
    secondary_attribute : str
        Name of a secondary instance attribute which supplies the value if
        ``default is None`` and no value was supplied during initialization.
    """

    def __init__(self, default=None, allowed_classes=None,
                 secondary_attribute=''):

        if allowed_classes is not None:
            self.allowed_classes = tuple(allowed_classes)
        else:
            self.allowed_classes = BaseDifferential

        super(DifferentialAttribute, self).__init__(default,
                                                    secondary_attribute)

    def convert_input(self, value):
        """
        Checks that the input is a differential object and is one of the
        allowed class types.

        Parameters
        ----------
        value : object
            Input value.

        Returns
        -------
        out, converted : correctly-typed object, boolean
            Tuple consisting of the correctly-typed object and a boolean which
            indicates if conversion was actually performed.

        Raises
        ------
        ValueError
            If the input is not valid for this attribute.
        """

        if not isinstance(value, self.allowed_classes):
            raise TypeError('Tried to set a DifferentialAttribute with '
                            'an unsupported Differential type {0}. Allowed '
                            'classes are: {1}'
                            .format(value.__class__,
                                    self.allowed_classes))

        return value, True


# Backwards-compatibility: these are the only classes that were previously
# released in v1.3
class FrameAttribute(Attribute):

    def __init__(self, *args, **kwargs):
        warnings.warn("FrameAttribute has been renamed to Attribute.",
                      AstropyDeprecationWarning)
        super(FrameAttribute, self).__init__(*args, **kwargs)

class TimeFrameAttribute(TimeAttribute):

    def __init__(self, *args, **kwargs):
        warnings.warn("TimeFrameAttribute has been renamed to TimeAttribute.",
                      AstropyDeprecationWarning)
        super(TimeFrameAttribute, self).__init__(*args, **kwargs)

class QuantityFrameAttribute(QuantityAttribute):

    def __init__(self, *args, **kwargs):
        warnings.warn("QuantityFrameAttribute has been renamed to "
                      "QuantityAttribute.", AstropyDeprecationWarning)
        super(QuantityFrameAttribute, self).__init__(*args, **kwargs)

class CartesianRepresentationFrameAttribute(CartesianRepresentationAttribute):

    def __init__(self, *args, **kwargs):
        warnings.warn("CartesianRepresentationFrameAttribute has been renamed "
                      "to CartesianRepresentationAttribute.",
                      AstropyDeprecationWarning)
        super(CartesianRepresentationFrameAttribute, self).__init__(
            *args, **kwargs)


# do this here to prevent a series of complicated circular imports
from .earth import EarthLocation
from .representation import CartesianRepresentation, BaseDifferential
