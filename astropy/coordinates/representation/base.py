# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Base classes for representations and differentials."""

import abc
import functools
import operator
import warnings
from typing import ClassVar, Final

import numpy as np

import astropy.units as u
from astropy.coordinates.angles import Angle
from astropy.utils import classproperty
from astropy.utils.data_info import MixinInfo
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import DuplicateRepresentationWarning
from astropy.utils.masked import MaskableShapedLikeNDArray, Masked, combine_masks

# Module-level dict mapping representation string alias names to classes.
# This is populated by __init_subclass__ when called by Representation or
# Differential classes so that they are all registered automatically.
REPRESENTATION_CLASSES = {}
DIFFERENTIAL_CLASSES = {}
# set for tracking duplicates
DUPLICATE_REPRESENTATIONS = set()


def _fqn_class(cls):
    """Get the fully qualified name of a class."""
    return cls.__module__ + "." + cls.__qualname__


@functools.cache
def get_reprdiff_cls_hash():
    """
    Returns a hash value that should be invariable if the
    `REPRESENTATION_CLASSES` and `DIFFERENTIAL_CLASSES` dictionaries have not
    changed.
    """
    return hash(tuple(REPRESENTATION_CLASSES.items())) + hash(
        tuple(DIFFERENTIAL_CLASSES.items())
    )


class BaseRepresentationOrDifferentialInfo(MixinInfo):
    """
    Container for meta information like name, description, format.  This is
    required when the object is used as a mixin column within a table, but can
    be used as a general way to store meta information.
    """

    attrs_from_parent = {"unit"}  # Indicates unit is read-only
    _supports_indexing = False
    mask_val = np.ma.masked

    @staticmethod
    def default_format(val):
        # Create numpy dtype so that numpy formatting will work.
        components = val.components
        values = tuple(getattr(val, component).value for component in components)
        a = np.empty(
            getattr(val, "shape", ()),
            [(component, value.dtype) for component, value in zip(components, values)],
        )
        for component, value in zip(components, values):
            a[component] = value
        return str(a)

    @property
    def _represent_as_dict_attrs(self):
        return self._parent.components

    @property
    def unit(self):
        if self._parent is None:
            return None

        unit = self._parent._unitstr
        return unit[1:-1] if unit.startswith("(") else unit

    def new_like(self, reps, length, metadata_conflicts="warn", name=None):
        """
        Return a new instance like ``reps`` with ``length`` rows.

        This is intended for creating an empty column object whose elements can
        be set in-place for table operations like join or vstack.

        Parameters
        ----------
        reps : list
            List of input representations or differentials.
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : `~astropy.coordinates.BaseRepresentation` or `~astropy.coordinates.BaseDifferential` subclass instance
            Empty instance of this class consistent with ``cols``

        """
        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(
            reps, metadata_conflicts, name, ("meta", "description")
        )

        # Make a new representation or differential with the desired length.
        rep0 = reps[0]
        out = rep0._apply(np.zeros_like, shape=(length,) + rep0.shape[1:])

        # Use __setitem__ machinery to check whether all representations
        # can represent themselves as this one without loss of information.
        # We use :0 to ensure we do not break on empty coordinates (with the
        # side benefit that we do not actually set anything).
        for rep in reps[1:]:
            try:
                out[:0] = rep[:0]
            except Exception as err:
                raise ValueError("input representations are inconsistent.") from err

        # Set (merged) info attributes.
        for attr in ("name", "meta", "description"):
            if attr in attrs:
                setattr(out.info, attr, attrs[attr])

        return out


class BaseRepresentationOrDifferential(MaskableShapedLikeNDArray):
    """3D coordinate representations and differentials.

    Parameters
    ----------
    comp1, comp2, comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D point or differential.  The names are the
        keys and the subclasses the values of the ``attr_classes`` attribute.
    copy : bool, optional
        If `True` (default), arrays will be copied; if `False`, they will be
        broadcast together but not use new memory.
    """

    # Ensure multiplication/division with ndarray or Quantity doesn't lead to
    # object arrays.
    __array_priority__: Final = 50000

    # Have to define this default b/c init_subclass is not called for the base class.
    name: ClassVar[str] = "base"
    """Name of the representation or differential.

    When a subclass is defined, by default, the name is the lower-cased name of the
    class with with any trailing 'representation' or 'differential' removed. (E.g.,
    'spherical' for `~astropy.coordinates.SphericalRepresentation` or
    `~astropy.coordinates.SphericalDifferential`.)

    This can be customized when defining a subclass by setting the class attribute.
    """

    info = BaseRepresentationOrDifferentialInfo()

    def __init_subclass__(cls) -> None:
        # Name of the representation or differential
        if "name" not in cls.__dict__:
            cls.name = (
                cls.__name__.lower()
                .removesuffix("representation")
                .removesuffix("differential")
            )

    def __init__(self, *args, **kwargs):
        # make argument a list, so we can pop them off.
        args = list(args)
        components = self.components
        if (
            args
            and isinstance(args[0], self.__class__)
            and all(arg is None for arg in args[1:])
        ):
            rep_or_diff = args[0]
            copy = kwargs.pop("copy", True)
            attrs = [getattr(rep_or_diff, component) for component in components]
            if "info" in rep_or_diff.__dict__:
                self.info = rep_or_diff.info

            if kwargs:
                raise TypeError(
                    "unexpected keyword arguments for case "
                    f"where class instance is passed in: {kwargs}"
                )

        else:
            attrs = []
            for component in components:
                try:
                    attr = args.pop(0) if args else kwargs.pop(component)
                except KeyError:
                    raise TypeError(
                        "__init__() missing 1 required positional "
                        f"argument: {component!r}"
                    ) from None

                if attr is None:
                    raise TypeError(
                        "__init__() missing 1 required positional argument:"
                        f" {component!r} (or first argument should be an instance of"
                        f" {self.__class__.__name__})."
                    )

                attrs.append(attr)

            copy = args.pop(0) if args else kwargs.pop("copy", True)

            if args:
                raise TypeError(f"unexpected arguments: {args}")

            if kwargs:
                for component in components:
                    if component in kwargs:
                        raise TypeError(
                            f"__init__() got multiple values for argument {component!r}"
                        )

                raise TypeError(f"unexpected keyword arguments: {kwargs}")

        # Pass attributes through the required initializing classes.
        attrs = [
            self.attr_classes[component](attr, copy=copy, subok=True)
            for component, attr in zip(components, attrs)
        ]
        try:
            bc_attrs = np.broadcast_arrays(*attrs, subok=True)
        except ValueError as err:
            if len(components) <= 2:
                c_str = " and ".join(components)
            else:
                c_str = ", ".join(components[:2]) + ", and " + components[2]
            raise ValueError(f"Input parameters {c_str} cannot be broadcast") from err

        # The output of np.broadcast_arrays() has limitations on writeability, so we perform
        # additional handling to enable writeability in most situations.  This is primarily
        # relevant for allowing the changing of the wrap angle of longitude components.
        #
        # If the shape has changed for a given component, broadcasting is needed:
        #     If copy=True, we make a copy of the broadcasted array to ensure writeability.
        #         Note that array had already been copied prior to the broadcasting.
        #         TODO: Find a way to avoid the double copy.
        #     If copy=False, we use the broadcasted array, and writeability may still be
        #         limited.
        # If the shape has not changed for a given component, we can proceed with using the
        #     non-broadcasted array, which avoids writeability issues from np.broadcast_arrays().
        attrs = [
            (bc_attr.copy() if copy else bc_attr)
            if bc_attr.shape != attr.shape
            else attr
            for attr, bc_attr in zip(attrs, bc_attrs)
        ]

        # Set private attributes for the attributes. (If not defined explicitly
        # on the class, the metaclass will define properties to access these.)
        for component, attr in zip(components, attrs):
            setattr(self, "_" + component, attr)

        # If any attribute has a mask, ensure all attributes are Masked.
        if any(hasattr(attr, "mask") for attr in attrs):
            self._ensure_masked()

    @deprecated("v7.1", alternative="name")
    @classmethod
    def get_name(cls):
        """Name of the representation or differential.

        Returns the ``.name`` attribute.
        """
        return cls.name

    # The two methods that any subclass has to define.
    @classmethod
    @abc.abstractmethod
    def from_cartesian(cls, other):
        """Create a representation of this class from a supplied Cartesian one.

        Parameters
        ----------
        other : `~astropy.coordinates.CartesianRepresentation`
            The representation to turn into this class

        Returns
        -------
        representation : `~astropy.coordinates.BaseRepresentation` subclass instance
            A new representation of this class's type.
        """
        # Note: the above docstring gets overridden for differentials.
        raise NotImplementedError()

    @abc.abstractmethod
    def to_cartesian(self):
        """Convert the representation to its Cartesian form.

        Note that any differentials get dropped.
        Also note that orientation information at the origin is *not* preserved by
        conversions through Cartesian coordinates. For example, transforming
        an angular position defined at distance=0 through cartesian coordinates
        and back will lose the original angular coordinates::

            >>> import astropy.units as u
            >>> import astropy.coordinates as coord
            >>> rep = coord.SphericalRepresentation(
            ...     lon=15*u.deg,
            ...     lat=-11*u.deg,
            ...     distance=0*u.pc)
            >>> rep.to_cartesian().represent_as(coord.SphericalRepresentation)
            <SphericalRepresentation (lon, lat, distance) in (rad, rad, pc)
                (0., 0., 0.)>

        Returns
        -------
        cartrepr : `~astropy.coordinates.CartesianRepresentation`
            The representation in Cartesian form.
        """
        # Note: the above docstring gets overridden for differentials.
        raise NotImplementedError()

    @property
    def components(self):
        """A tuple with the in-order names of the coordinate components."""
        return tuple(self.attr_classes)

    def __eq__(self, value):
        """Equality operator.

        This implements strict equality and requires that the representation
        classes are identical and that the representation data are exactly equal.
        """
        if self.__class__ is not value.__class__:
            raise TypeError(
                "cannot compare: objects must have same class: "
                f"{self.__class__.__name__} vs. {value.__class__.__name__}"
            )

        try:
            np.broadcast(self, value)
        except ValueError as exc:
            raise ValueError(f"cannot compare: {exc}") from exc

        out = True
        for comp in self.components:
            out &= getattr(self, "_" + comp) == getattr(value, "_" + comp)

        return out

    def __ne__(self, value):
        return np.logical_not(self == value)

    def _apply(self, method, *args, **kwargs):
        """Create a new representation or differential with ``method`` applied
        to the component data.

        In typical usage, the method is any of the shape-changing methods for
        `~numpy.ndarray` (``reshape``, ``swapaxes``, etc.), as well as those
        picking particular elements (``__getitem__``, ``take``, etc.), which
        are all defined in `~astropy.utils.shapes.ShapedLikeNDArray`. It will be
        applied to the underlying arrays (e.g., ``x``, ``y``, and ``z`` for
        `~astropy.coordinates.CartesianRepresentation`), with the results used
        to create a new instance.

        Internally, it is also used to apply functions to the components
        (in particular, `~numpy.broadcast_to`).

        Parameters
        ----------
        method : str or callable
            If str, it is the name of a method that is applied to the internal
            ``components``. If callable, the function is applied.
        *args : tuple
            Any positional arguments for ``method``.
        **kwargs : dict
            Any keyword arguments for ``method``.
        """
        if callable(method):
            apply_method = lambda array: method(array, *args, **kwargs)
        else:
            apply_method = operator.methodcaller(method, *args, **kwargs)

        new = super().__new__(self.__class__)
        for component in self.components:
            setattr(new, "_" + component, apply_method(getattr(self, component)))

        # Copy other 'info' attr only if it has actually been defined.
        # See PR #3898 for further explanation and justification, along
        # with Quantity.__array_finalize__
        if "info" in self.__dict__:
            new.info = self.info

        return new

    def __setitem__(self, item, value):
        set_mask = value is np.ma.masked
        clear_mask = value is np.ma.nomask
        if not (value.__class__ is self.__class__ or set_mask or clear_mask):
            raise TypeError(
                "can only set from object of same class: "
                f"{self.__class__.__name__} vs. {value.__class__.__name__}"
                " (unless setting or clearing the mask with"
                " np.ma.masked or np.ma.nomask)."
            )

        if not self.masked:
            if clear_mask:
                # Clearing masked elements on an unmasked instance: nothing to do.
                return

            # Ensure our components are masked if a mask needs to be set.
            # NOTE: we could also make ourselves masked if value.masked.
            # But then we have to be sure that Time does the same, and live
            # with the inconsistency that things like ndarray and Quantity cannot
            # become masked when setting an item with a masked value.  See
            # https://github.com/astropy/astropy/pull/17016#issuecomment-2439607869
            if set_mask:
                self._ensure_masked()

        if set_mask or clear_mask:
            for comp in self.components:
                c = "_" + comp
                getattr(self, c).mask[item] = set_mask
            return

        for component in self.components:
            c = "_" + component
            getattr(self, c)[item] = getattr(value, c)

    @property
    def shape(self):
        """The shape of the instance and underlying arrays.

        Like `~numpy.ndarray.shape`, can be set to a new shape by assigning a
        tuple.  Note that if different instances share some but not all
        underlying data, setting the shape of one instance can make the other
        instance unusable.  Hence, it is strongly recommended to get new,
        reshaped instances with the ``reshape`` method.

        Raises
        ------
        ValueError
            If the new shape has the wrong total number of elements.
        AttributeError
            If the shape of any of the components cannot be changed without the
            arrays being copied.  For these cases, use the ``reshape`` method
            (which copies any arrays that cannot be reshaped in-place).
        """
        return getattr(self, self.components[0]).shape

    @shape.setter
    def shape(self, shape):
        # We keep track of arrays that were already reshaped since we may have
        # to return those to their original shape if a later shape-setting
        # fails. (This can happen since coordinates are broadcast together.)
        reshaped = []
        oldshape = self.shape
        for component in self.components:
            val = getattr(self, component)
            if val.size > 1:
                try:
                    val.shape = shape
                except Exception:
                    for val2 in reshaped:
                        val2.shape = oldshape
                    raise
                else:
                    reshaped.append(val)

    @property
    def masked(self):
        return isinstance(getattr(self, self.components[0]), Masked)

    def _ensure_masked(self):
        """Ensure Masked components."""
        # TODO: should we just allow the above property to be set?
        # But be sure the API remains consistent with Time!
        if not self.masked:
            for comp in self.components:
                c = "_" + comp
                setattr(self, c, Masked(getattr(self, c)))

    def get_mask(self, *attrs):
        """Calculate the mask, by combining masks from the given attributes.

        Parameters
        ----------
        *attrs : str
            Attributes from which to get the masks to combine. If not given,
            use all components of the class.

        Returns
        -------
        mask : ~numpy.ndarray of bool
            The combined, read-only mask. If the instance is not masked, it
            is an array of `False` with the correct shape.
        """
        if not attrs:
            attrs = self.components

        values = operator.attrgetter(*attrs)(self)
        if not isinstance(values, tuple):
            values = (values,)

        mask = combine_masks([getattr(v, "mask", None) for v in values])
        return np.broadcast_to(mask, self.shape)  # Makes it readonly too.

    mask = property(get_mask, doc="The combined mask of all components.")

    # Required to support multiplication and division, and defined by the base
    # representation and differential classes.
    @abc.abstractmethod
    def _scale_operation(self, op, *args):
        raise NotImplementedError()

    def __mul__(self, other):
        return self._scale_operation(operator.mul, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._scale_operation(operator.truediv, other)

    def __neg__(self):
        return self._scale_operation(operator.neg)

    # Follow numpy convention and make an independent copy.
    def __pos__(self):
        return self.copy()

    # Required to support addition and subtraction, and defined by the base
    # representation and differential classes.
    @abc.abstractmethod
    def _combine_operation(self, op, other, reverse=False):
        raise NotImplementedError()

    def __add__(self, other):
        return self._combine_operation(operator.add, other)

    def __radd__(self, other):
        return self._combine_operation(operator.add, other, reverse=True)

    def __sub__(self, other):
        return self._combine_operation(operator.sub, other)

    def __rsub__(self, other):
        return self._combine_operation(operator.sub, other, reverse=True)

    # The following are used for repr and str
    @property
    def _values(self):
        """Turn the coordinates into a record array with the coordinate values.

        The record array fields will have the component names.
        """
        coo_items = [(c, getattr(self, c)) for c in self.components]
        result = np.empty_like(
            coo_items[0][1].value, dtype=[(c, coo.dtype) for c, coo in coo_items]
        )
        for c, coo in coo_items:
            result[c] = coo.value
        return result

    @property
    def _units(self):
        """Return a dictionary with the units of the coordinate components."""
        return {cmpnt: getattr(self, cmpnt).unit for cmpnt in self.components}

    @property
    def _unitstr(self):
        units = self._units.values()
        if len(units_set := set(units)) == 1:
            return str(units_set.pop())
        return f"({', '.join(map(str, units))})"

    def __str__(self):
        return f"{np.array2string(self._values, separator=', ')} {self._unitstr:s}"

    def __repr__(self):
        prefixstr = "    "
        arrstr = np.array2string(self._values, prefix=prefixstr, separator=", ")

        diffstr = ""
        if diffs := getattr(self, "differentials", None):
            diffstr = f"\n (has differentials w.r.t.: {', '.join(map(repr, diffs))})"

        unitstr = ("in " + self._unitstr) if self._unitstr else "[dimensionless]"
        return (
            f"<{self.__class__.__name__} ({', '.join(self.components)})"
            f" {unitstr:s}\n{prefixstr}{arrstr}{diffstr}>"
        )


class RepresentationInfo(BaseRepresentationOrDifferentialInfo):
    @property
    def _represent_as_dict_attrs(self):
        attrs = super()._represent_as_dict_attrs
        if self._parent._differentials:
            attrs += ("differentials",)
        return attrs

    def _represent_as_dict(self, attrs=None):
        out = super()._represent_as_dict(attrs)
        for key, value in out.pop("differentials", {}).items():
            out[f"differentials.{key}"] = value
        return out

    def _construct_from_dict(self, map):
        differentials = {}
        for key in list(map.keys()):
            if key.startswith("differentials."):
                differentials[key[14:]] = map.pop(key)
        map["differentials"] = differentials
        return super()._construct_from_dict(map)


class BaseRepresentation(BaseRepresentationOrDifferential):
    """Base for representing a point in a 3D coordinate system.

    Parameters
    ----------
    comp1, comp2, comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D points.  The names are the keys and the
        subclasses the values of the ``attr_classes`` attribute.
    differentials : dict, `~astropy.coordinates.BaseDifferential`, optional
        Any differential classes that should be associated with this
        representation. The input must either be a single `~astropy.coordinates.BaseDifferential`
        subclass instance, or a dictionary with keys set to a string
        representation of the SI unit with which the differential (derivative)
        is taken. For example, for a velocity differential on a positional
        representation, the key would be ``'s'`` for seconds, indicating that
        the derivative is a time derivative.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.

    Notes
    -----
    All representation classes should subclass this base representation class,
    and define an ``attr_classes`` attribute, a `dict`
    which maps component names to the class that creates them. They must also
    define a ``to_cartesian`` method and a ``from_cartesian`` class method. By
    default, transformations are done via the cartesian system, but classes
    that want to define a smarter transformation path can overload the
    ``represent_as`` method. If one wants to use an associated differential
    class, one should also define ``unit_vectors`` and ``scale_factors``
    methods (see those methods for details).
    """

    info = RepresentationInfo()
    # Ensure _differentials always exists.
    _differentials = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)  # sets `cls.name`

        # Register representation name (except for bases on which other
        # representations are built, but which cannot themselves be used).
        if cls.__name__.startswith("Base"):
            return

        if not hasattr(cls, "attr_classes"):
            raise NotImplementedError(
                'Representations must have an "attr_classes" class attribute.'
            )

        repr_name = cls.name
        # first time a duplicate is added
        # remove first entry and add both using their qualnames
        if repr_name in REPRESENTATION_CLASSES:
            DUPLICATE_REPRESENTATIONS.add(repr_name)

            fqn_cls = _fqn_class(cls)
            existing = REPRESENTATION_CLASSES[repr_name]
            fqn_existing = _fqn_class(existing)

            if fqn_cls == fqn_existing:
                raise ValueError(f'Representation "{fqn_cls}" already defined')

            msg = (
                f'Representation "{repr_name}" already defined, removing it to avoid'
                f' confusion.Use qualnames "{fqn_cls}" and "{fqn_existing}" or class'
                " instances directly"
            )
            warnings.warn(msg, DuplicateRepresentationWarning)

            del REPRESENTATION_CLASSES[repr_name]
            REPRESENTATION_CLASSES[fqn_existing] = existing
            repr_name = fqn_cls

        # further definitions with the same name, just add qualname
        elif repr_name in DUPLICATE_REPRESENTATIONS:
            fqn_cls = _fqn_class(cls)
            warnings.warn(
                f'Representation "{repr_name}" already defined, using qualname '
                f'"{fqn_cls}".'
            )
            repr_name = fqn_cls
            if repr_name in REPRESENTATION_CLASSES:
                raise ValueError(f'Representation "{repr_name}" already defined')

        REPRESENTATION_CLASSES[repr_name] = cls
        get_reprdiff_cls_hash.cache_clear()

        # define getters for any component that does not yet have one.
        for component in cls.attr_classes:
            if not hasattr(cls, component):
                setattr(
                    cls,
                    component,
                    property(
                        lambda self, comp=f"_{component}": getattr(self, comp),
                        doc=f"The '{component}' component of the points(s).",
                    ),
                )

    def __init__(self, *args, differentials=None, **kwargs):
        # Handle any differentials passed in.
        super().__init__(*args, **kwargs)
        if differentials is None and args and isinstance(args[0], self.__class__):
            differentials = args[0]._differentials
        self._differentials = self._validate_differentials(differentials)
        # If any part is masked, all should be.
        if self.masked or any(d.masked for d in self._differentials.values()):
            self._ensure_masked()

    def _ensure_masked(self):
        super()._ensure_masked()
        for d in self._differentials.values():
            d._ensure_masked()

    def _validate_differentials(self, differentials):
        """
        Validate that the provided differentials are appropriate for this
        representation and recast/reshape as necessary and then return.

        Note that this does *not* set the differentials on
        ``self._differentials``, but rather leaves that for the caller.
        """
        from .spherical import RadialDifferential, UnitSphericalRepresentation

        # Now handle the actual validation of any specified differential classes
        if differentials is None:
            differentials = {}

        elif isinstance(differentials, BaseDifferential):
            # We can't handle auto-determining the key for this combo
            if isinstance(differentials, RadialDifferential) and isinstance(
                self, UnitSphericalRepresentation
            ):
                raise ValueError(
                    "To attach a RadialDifferential to a UnitSphericalRepresentation,"
                    " you must supply a dictionary with an appropriate key."
                )

            key = differentials._get_deriv_key(self)
            differentials = {key: differentials}

        for key in differentials:
            try:
                diff = differentials[key]
            except TypeError as err:
                raise TypeError(
                    "'differentials' argument must be a dictionary-like object"
                ) from err

            diff._check_base(self)

            if isinstance(diff, RadialDifferential) and isinstance(
                self, UnitSphericalRepresentation
            ):
                # We trust the passing of a key for a RadialDifferential
                # attached to a UnitSphericalRepresentation because it will not
                # have a paired component name (UnitSphericalRepresentation has
                # no .distance) to automatically determine the expected key
                pass

            else:
                expected_key = diff._get_deriv_key(self)
                if key != expected_key:
                    raise ValueError(
                        f"For differential object '{repr(diff)}', expected "
                        f"unit key = '{expected_key}' but received key = '{key}'"
                    )

            # For now, we are very rigid: differentials must have the same shape
            # as the representation. This makes it easier to handle __getitem__
            # and any other shape-changing operations on representations that
            # have associated differentials
            if diff.shape != self.shape:
                # TODO: message of IncompatibleShapeError is not customizable,
                #       so use a valueerror instead?
                raise ValueError(
                    "Shape of differentials must be the same "
                    f"as the shape of the representation ({diff.shape} vs {self.shape})"
                )

        return differentials

    def _raise_if_has_differentials(self, op_name):
        """
        Used to raise a consistent exception for any operation that is not
        supported when a representation has differentials attached.
        """
        if self.differentials:
            raise TypeError(
                f"Operation '{op_name}' is not supported when "
                f"differentials are attached to a {self.__class__.__name__}."
            )

    @classproperty
    def _compatible_differentials(cls):
        return [DIFFERENTIAL_CLASSES[cls.name]]

    @property
    def differentials(self):
        """A dictionary of differential class instances.

        The keys of this dictionary must be a string representation of the SI
        unit with which the differential (derivative) is taken. For example, for
        a velocity differential on a positional representation, the key would be
        ``'s'`` for seconds, indicating that the derivative is a time
        derivative.
        """
        return self._differentials

    # We do not make unit_vectors and scale_factors abstract methods, since
    # they are only necessary if one also defines an associated Differential.
    # Also, doing so would break pre-differential representation subclasses.
    def unit_vectors(self):
        r"""Cartesian unit vectors in the direction of each component.

        Given unit vectors :math:`\hat{e}_c` and scale factors :math:`f_c`,
        a change in one component of :math:`\delta c` corresponds to a change
        in representation of :math:`\delta c \times f_c \times \hat{e}_c`.

        Returns
        -------
        unit_vectors : dict of `~astropy.coordinates.CartesianRepresentation`
            The keys are the component names.
        """
        raise NotImplementedError(f"{type(self)} has not implemented unit vectors")

    def scale_factors(self):
        r"""Scale factors for each component's direction.

        Given unit vectors :math:`\hat{e}_c` and scale factors :math:`f_c`,
        a change in one component of :math:`\delta c` corresponds to a change
        in representation of :math:`\delta c \times f_c \times \hat{e}_c`.

        Returns
        -------
        scale_factors : dict of `~astropy.units.Quantity`
            The keys are the component names.
        """
        raise NotImplementedError(f"{type(self)} has not implemented scale factors.")

    def _re_represent_differentials(self, new_rep, differential_class):
        """Re-represent the differentials to the specified classes.

        This returns a new dictionary with the same keys but with the
        attached differentials converted to the new differential classes.
        """
        if differential_class is None:
            return {}

        if not self.differentials and differential_class:
            raise ValueError("No differentials associated with this representation!")

        if (
            len(self.differentials) == 1
            and isinstance(differential_class, type)
            and issubclass(differential_class, BaseDifferential)
        ):
            differential_class = {
                next(iter(self.differentials.keys())): differential_class
            }

        elif differential_class.keys() != self.differentials.keys():
            raise ValueError(
                "Desired differential classes must be passed in as a dictionary with"
                " keys equal to a string representation of the unit of the derivative"
                " for each differential stored with this "
                f"representation object ({self.differentials})"
            )

        new_diffs = {}
        for k in self.differentials:
            diff = self.differentials[k]
            try:
                new_diffs[k] = diff.represent_as(differential_class[k], base=self)
            except Exception as err:
                if differential_class[k] not in new_rep._compatible_differentials:
                    raise TypeError(
                        f"Desired differential class {differential_class[k]} is not "
                        "compatible with the desired "
                        f"representation class {new_rep.__class__}"
                    ) from err
                raise

        return new_diffs

    def represent_as(self, other_class, differential_class=None):
        """Convert coordinates to another representation.

        If the instance is of the requested class, it is returned unmodified.
        By default, conversion is done via Cartesian coordinates.
        Also note that orientation information at the origin is *not* preserved by
        conversions through Cartesian coordinates. See the docstring for
        :meth:`~astropy.coordinates.BaseRepresentationOrDifferential.to_cartesian`
        for an example.

        Parameters
        ----------
        other_class : `~astropy.coordinates.BaseRepresentation` subclass
            The type of representation to turn the coordinates into.
        differential_class : dict of `~astropy.coordinates.BaseDifferential`, optional
            Classes in which the differentials should be represented.
            Can be a single class if only a single differential is attached,
            otherwise it should be a `dict` keyed by the same keys as the
            differentials.
        """
        if other_class is self.__class__ and not differential_class:
            return self.without_differentials()

        else:
            if isinstance(other_class, str):
                raise ValueError(
                    "Input to a representation's represent_as must be a class, not "
                    "a string. For strings, use frame objects."
                )

            if other_class is not self.__class__:
                # The default is to convert via cartesian coordinates
                new_rep = other_class.from_cartesian(self.to_cartesian())
            else:
                new_rep = self

            new_rep._differentials = self._re_represent_differentials(
                new_rep, differential_class
            )

            return new_rep

    def transform(self, matrix):
        """Transform coordinates using a 3x3 matrix in a Cartesian basis.

        This returns a new representation and does not modify the original one.
        Any differentials attached to this representation will also be
        transformed.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 (or stack thereof) matrix, such as a rotation matrix.

        """
        from .cartesian import CartesianDifferential, CartesianRepresentation

        # route transformation through Cartesian
        crep = self.represent_as(
            CartesianRepresentation,
            differential_class=dict.fromkeys(self.differentials, CartesianDifferential),
        ).transform(matrix)

        # move back to original representation
        difs_cls = {k: diff.__class__ for k, diff in self.differentials.items()}
        rep = crep.represent_as(self.__class__, difs_cls)
        return rep

    def with_differentials(self, differentials):
        """
        Create a new representation with the same positions as this
        representation, but with these new differentials.

        Differential keys that already exist in this object's differential dict
        are overwritten.

        Parameters
        ----------
        differentials : sequence of `~astropy.coordinates.BaseDifferential` subclass instance
            The differentials for the new representation to have.

        Returns
        -------
        `~astropy.coordinates.BaseRepresentation` subclass instance
            A copy of this representation, but with the ``differentials`` as
            its differentials.
        """
        if not differentials:
            return self

        differentials = self._validate_differentials(differentials)
        return self.__class__(
            *[getattr(self, component) for component in self.components],
            differentials=self.differentials | differentials,
            copy=False,
        )

    def without_differentials(self):
        """Return a copy of the representation without attached differentials.

        Returns
        -------
        `~astropy.coordinates.BaseRepresentation` subclass instance
            A shallow copy of this representation, without any differentials.
            If no differentials were present, no copy is made.
        """
        if not self._differentials:
            return self

        args = [getattr(self, component) for component in self.components]
        return self.__class__(*args, copy=False)

    @classmethod
    def from_representation(cls, representation):
        """Create a new instance of this representation from another one.

        Parameters
        ----------
        representation : `~astropy.coordinates.BaseRepresentation` instance
            The presentation that should be converted to this class.
        """
        return representation.represent_as(cls)

    def __eq__(self, value):
        """Equality operator for BaseRepresentation.

        This implements strict equality and requires that the representation
        classes are identical, the differentials are identical, and that the
        representation data are exactly equal.
        """
        # BaseRepresentationOrDifferental (checks classes and compares components)
        out = super().__eq__(value)

        # super() checks that the class is identical so can this even happen?
        # (same class, different differentials ?)
        if self._differentials.keys() != value._differentials.keys():
            raise ValueError("cannot compare: objects must have same differentials")

        for self_diff, value_diff in zip(
            self._differentials.values(), value._differentials.values()
        ):
            out &= self_diff == value_diff

        return out

    def __ne__(self, value):
        return np.logical_not(self == value)

    def _apply(self, method, *args, **kwargs):
        """Create a new representation with ``method`` applied to the component
        data.

        This is not a simple inherit from ``BaseRepresentationOrDifferential``
        because we need to call ``._apply()`` on any associated differential
        classes.

        See docstring for `BaseRepresentationOrDifferential._apply`.

        Parameters
        ----------
        method : str or callable
            If str, it is the name of a method that is applied to the internal
            ``components``. If callable, the function is applied.
        *args : tuple
            Any positional arguments for ``method``.
        **kwargs : dict
            Any keyword arguments for ``method``.

        """
        rep = super()._apply(method, *args, **kwargs)

        rep._differentials = {
            k: diff._apply(method, *args, **kwargs)
            for k, diff in self._differentials.items()
        }
        return rep

    def __setitem__(self, item, value):
        if value is np.ma.masked or value is np.ma.nomask:
            return super().__setitem__(item, value)

        if not isinstance(value, BaseRepresentation):
            raise TypeError(
                f"value must be a representation instance, not {type(value)}."
            )

        if not (
            isinstance(value, self.__class__)
            or len(value.attr_classes) == len(self.attr_classes)
        ):
            raise ValueError(
                f"value must be representable as {self.__class__.__name__} "
                "without loss of information."
            )

        diff_classes = {}
        if self._differentials:
            if self._differentials.keys() != value._differentials.keys():
                raise ValueError("value must have the same differentials.")

            for key, self_diff in self._differentials.items():
                diff_classes[key] = self_diff_cls = self_diff.__class__
                value_diff_cls = value._differentials[key].__class__
                if not (
                    isinstance(value_diff_cls, self_diff_cls)
                    or (
                        len(value_diff_cls.attr_classes)
                        == len(self_diff_cls.attr_classes)
                    )
                ):
                    raise ValueError(
                        f"value differential {key!r} must be representable as "
                        f"{self_diff.__class__.__name__} without loss of information."
                    )

        value = value.represent_as(self.__class__, diff_classes)
        super().__setitem__(item, value)
        for key, differential in self._differentials.items():
            differential[item] = value._differentials[key]

    def _scale_operation(self, op, *args):
        """Scale all non-angular components, leaving angular ones unchanged.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.mul`, `~operator.neg`, etc.
        *args
            Any arguments required for the operator (typically, what is to
            be multiplied with, divided by).
        """
        results = []
        for component, cls in self.attr_classes.items():
            value = getattr(self, component)
            if issubclass(cls, Angle):
                results.append(value)
            else:
                results.append(op(value, *args))

        # try/except catches anything that cannot initialize the class, such
        # as operations that returned NotImplemented or a representation
        # instead of a quantity (as would happen for, e.g., rep * rep).
        try:
            result = self.__class__(*results)
        except Exception:
            return NotImplemented

        for key, differential in self.differentials.items():
            diff_result = differential._scale_operation(op, *args, scaled_base=True)
            result.differentials[key] = diff_result

        return result

    def _combine_operation(self, op, other, reverse=False):
        """Combine two representation.

        By default, operate on the cartesian representations of both.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            The other representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        self._raise_if_has_differentials(op.__name__)

        result = self.to_cartesian()._combine_operation(op, other, reverse)
        if result is NotImplemented:
            return NotImplemented
        else:
            return self.from_cartesian(result)

    # We need to override this setter to support differentials
    @BaseRepresentationOrDifferential.shape.setter
    def shape(self, shape):
        orig_shape = self.shape

        # See: https://stackoverflow.com/questions/3336767/ for an example
        BaseRepresentationOrDifferential.shape.fset(self, shape)

        # also try to perform shape-setting on any associated differentials
        try:
            for k in self.differentials:
                self.differentials[k].shape = shape
        except Exception:
            BaseRepresentationOrDifferential.shape.fset(self, orig_shape)
            for k in self.differentials:
                self.differentials[k].shape = orig_shape

            raise

    def norm(self):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Note that any associated differentials will be dropped during this
        operation.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        return np.sqrt(
            sum(
                getattr(self, component) ** 2
                for component, cls in self.attr_classes.items()
                if not issubclass(cls, Angle)
            )
        )

    def mean(self, *args, **kwargs):
        """Vector mean.

        Averaging is done by converting the representation to cartesian, and
        taking the mean of the x, y, and z components. The result is converted
        back to the same representation as the input.

        Refer to `~numpy.mean` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.

        Returns
        -------
        mean : `~astropy.coordinates.BaseRepresentation` subclass instance
            Vector mean, in the same representation as that of the input.
        """
        self._raise_if_has_differentials("mean")
        return self.from_cartesian(self.to_cartesian().mean(*args, **kwargs))

    def sum(self, *args, **kwargs):
        """Vector sum.

        Adding is done by converting the representation to cartesian, and
        summing the x, y, and z components. The result is converted back to the
        same representation as the input.

        Refer to `~numpy.sum` for full documentation of the arguments, noting
        that ``axis`` is the entry in the ``shape`` of the representation, and
        that the ``out`` argument cannot be used.

        Returns
        -------
        sum : `~astropy.coordinates.BaseRepresentation` subclass instance
            Vector sum, in the same representation as that of the input.
        """
        self._raise_if_has_differentials("sum")
        return self.from_cartesian(self.to_cartesian().sum(*args, **kwargs))

    def dot(self, other):
        """Dot product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`.

        Note that any associated differentials will be dropped during this
        operation.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseRepresentation`
            The representation to take the dot product with.

        Returns
        -------
        dot_product : `~astropy.units.Quantity`
            The sum of the product of the x, y, and z components of the
            cartesian representations of ``self`` and ``other``.
        """
        return self.to_cartesian().dot(other)

    def cross(self, other):
        """Vector cross product of two representations.

        The calculation is done by converting both ``self`` and ``other``
        to `~astropy.coordinates.CartesianRepresentation`, and converting the
        result back to the type of representation of ``self``.

        Parameters
        ----------
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            The representation to take the cross product with.

        Returns
        -------
        cross_product : `~astropy.coordinates.BaseRepresentation` subclass instance
            With vectors perpendicular to both ``self`` and ``other``, in the
            same type of representation as ``self``.
        """
        self._raise_if_has_differentials("cross")
        return self.from_cartesian(self.to_cartesian().cross(other))


class BaseDifferential(BaseRepresentationOrDifferential):
    r"""A base class representing differentials of representations.

    These represent differences or derivatives along each component.
    E.g., for physics spherical coordinates, these would be
    :math:`\delta r, \delta \theta, \delta \phi`.

    Parameters
    ----------
    d_comp1, d_comp2, d_comp3 : `~astropy.units.Quantity` or subclass
        The components of the 3D differentials.  The names are the keys and the
        subclasses the values of the ``attr_classes`` attribute.
    copy : bool, optional
        If `True` (default), arrays will be copied. If `False`, arrays will
        be references, though possibly broadcast to ensure matching shapes.

    Notes
    -----
    All differential representation classes should subclass this base class,
    and define an ``base_representation`` attribute with the class of the
    regular `~astropy.coordinates.BaseRepresentation` for which differential
    coordinates are provided. This will set up a default ``attr_classes``
    instance with names equal to the base component names prefixed by ``d_``,
    and all classes set to `~astropy.units.Quantity`, plus properties to access
    those, and a default ``__init__`` for initialization.
    """

    def __init_subclass__(cls, **kwargs):
        """Set default ``attr_classes`` and component getters on a Differential.

        For these, the components are those of the base representation prefixed
        by 'd_', and the class is `~astropy.units.Quantity`.
        """
        super().__init_subclass__(**kwargs)  # sets `cls.name`

        # Don't do anything for base helper classes.
        if cls.__name__ in (
            "BaseDifferential",
            "BaseSphericalDifferential",
            "BaseSphericalCosLatDifferential",
        ):
            return

        if not hasattr(cls, "base_representation"):
            raise NotImplementedError(
                "Differential representations must have a"
                '"base_representation" class attribute.'
            )

        # If not defined explicitly, create attr_classes.
        if not hasattr(cls, "attr_classes"):
            base_attr_classes = cls.base_representation.attr_classes
            cls.attr_classes = {"d_" + c: u.Quantity for c in base_attr_classes}

        repr_name = cls.name
        if repr_name in DIFFERENTIAL_CLASSES:
            raise ValueError(f"Differential class {repr_name} already defined")

        DIFFERENTIAL_CLASSES[repr_name] = cls
        get_reprdiff_cls_hash.cache_clear()

        # If not defined explicitly, create properties for the components.
        for component in cls.attr_classes:
            if not hasattr(cls, component):
                setattr(
                    cls,
                    component,
                    property(
                        lambda self, comp=f"_{component}": getattr(self, comp),
                        doc=f"Component '{component}' of the Differential.",
                    ),
                )

    @classmethod
    def _check_base(cls, base):
        if cls not in base._compatible_differentials:
            raise TypeError(
                f"Differential class {cls} is not compatible with the "
                f"base (representation) class {base.__class__}"
            )

    def _get_deriv_key(self, base):
        """Given a base (representation instance), determine the unit of the
        derivative by removing the representation unit from the component units
        of this differential.
        """
        # This check is just a last resort so we don't return a strange unit key
        # from accidentally passing in the wrong base.
        self._check_base(base)

        for name in base.components:
            comp = getattr(base, name)
            d_comp = getattr(self, f"d_{name}", None)
            if d_comp is not None:
                d_unit = comp.unit / d_comp.unit

                # This is quite a bit faster than using to_system() or going
                # through Quantity()
                d_unit_si = d_unit.decompose(u.si.bases)
                d_unit_si._scale = 1  # remove the scale from the unit

                return str(d_unit_si)
        raise RuntimeError(
            "Invalid representation-differential units! This likely happened "
            "because either the representation or the associated differential "
            "have non-standard units. Check that the input positional data have "
            "positional units, and the input velocity data have velocity units, "
            "or are both dimensionless."
        )

    @classmethod
    def _get_base_vectors(cls, base):
        """Get unit vectors and scale factors from base.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            The points for which the unit vectors and scale factors should be
            retrieved.

        Returns
        -------
        unit_vectors : dict of `~astropy.coordinates.CartesianRepresentation`
            In the directions of the coordinates of base.
        scale_factors : dict of `~astropy.units.Quantity`
            Scale factors for each of the coordinates

        Raises
        ------
        TypeError : if the base is not of the correct type
        """
        cls._check_base(base)
        return base.unit_vectors(), base.scale_factors()

    def to_cartesian(self, base):
        """Convert the differential to 3D rectangular cartesian coordinates.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            The points for which the differentials are to be converted: each of
            the components is multiplied by its unit vectors and scale factors.

        Returns
        -------
        `~astropy.coordinates.CartesianDifferential`
            This object, converted.

        """
        base_e, base_sf = self._get_base_vectors(base)
        return functools.reduce(
            operator.add,
            (
                getattr(self, d_c) * base_sf[c] * base_e[c]
                for d_c, c in zip(self.components, base.components)
            ),
        )

    @classmethod
    def from_cartesian(cls, other, base):
        """Convert the differential from 3D rectangular cartesian coordinates to
        the desired class.

        Parameters
        ----------
        other
            The object to convert into this differential.
        base : `~astropy.coordinates.BaseRepresentation`
            The points for which the differentials are to be converted: each of
            the components is multiplied by its unit vectors and scale factors.
            Will be converted to ``cls.base_representation`` if needed.

        Returns
        -------
        `~astropy.coordinates.BaseDifferential` subclass instance
            A new differential object that is this class' type.
        """
        base = base.represent_as(cls.base_representation)
        base_e, base_sf = cls._get_base_vectors(base)
        return cls(
            *(other.dot(e / base_sf[component]) for component, e in base_e.items()),
            copy=False,
        )

    def represent_as(self, other_class, base):
        """Convert coordinates to another representation.

        If the instance is of the requested class, it is returned unmodified.
        By default, conversion is done via cartesian coordinates.

        Parameters
        ----------
        other_class : `~astropy.coordinates.BaseRepresentation` subclass
            The type of representation to turn the coordinates into.
        base : instance of ``self.base_representation``
            Base relative to which the differentials are defined.  If the other
            class is a differential representation, the base will be converted
            to its ``base_representation``.
        """
        if other_class is self.__class__:
            return self

        # The default is to convert via cartesian coordinates.
        self_cartesian = self.to_cartesian(base)
        if issubclass(other_class, BaseDifferential):
            return other_class.from_cartesian(self_cartesian, base)
        else:
            return other_class.from_cartesian(self_cartesian)

    @classmethod
    def from_representation(cls, representation, base):
        """Create a new instance of this representation from another one.

        Parameters
        ----------
        representation : `~astropy.coordinates.BaseRepresentation` instance
            The presentation that should be converted to this class.
        base : instance of ``cls.base_representation``
            The base relative to which the differentials will be defined. If
            the representation is a differential itself, the base will be
            converted to its ``base_representation`` to help convert it.
        """
        if isinstance(representation, BaseDifferential):
            cartesian = representation.to_cartesian(
                base.represent_as(representation.base_representation)
            )
        else:
            cartesian = representation.to_cartesian()

        return cls.from_cartesian(cartesian, base)

    def transform(self, matrix, base, transformed_base):
        """Transform differential using a 3x3 matrix in a Cartesian basis.

        This returns a new differential and does not modify the original one.

        Parameters
        ----------
        matrix : (3,3) array-like
            A 3x3 (or stack thereof) matrix, such as a rotation matrix.
        base : instance of ``cls.base_representation``
            Base relative to which the differentials are defined.  If the other
            class is a differential representation, the base will be converted
            to its ``base_representation``.
        transformed_base : instance of ``cls.base_representation``
            Base relative to which the transformed differentials are defined.
            If the other class is a differential representation, the base will
            be converted to its ``base_representation``.
        """
        from .cartesian import CartesianDifferential

        # route transformation through Cartesian
        cdiff = self.represent_as(CartesianDifferential, base=base).transform(matrix)
        # move back to original representation
        diff = cdiff.represent_as(self.__class__, transformed_base)
        return diff

    def _scale_operation(self, op, *args, scaled_base=False):
        """Scale all components.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.mul`, `~operator.neg`, etc.
        *args
            Any arguments required for the operator (typically, what is to
            be multiplied with, divided by).
        scaled_base : bool, optional
            Whether the base was scaled the same way. This affects whether
            differential components should be scaled. For instance, a differential
            in longitude should not be scaled if its spherical base is scaled
            in radius.
        """
        scaled_attrs = [op(getattr(self, c), *args) for c in self.components]
        return self.__class__(*scaled_attrs, copy=False)

    def _combine_operation(self, op, other, reverse=False):
        """Combine two differentials, or a differential with a representation.

        If ``other`` is of the same differential type as ``self``, the
        components will simply be combined.  If ``other`` is a representation,
        it will be used as a base for which to evaluate the differential,
        and the result is a new representation.

        Parameters
        ----------
        op : `~operator` callable
            Operator to apply (e.g., `~operator.add`, `~operator.sub`, etc.
        other : `~astropy.coordinates.BaseRepresentation` subclass instance
            The other differential or representation.
        reverse : bool
            Whether the operands should be reversed (e.g., as we got here via
            ``self.__rsub__`` because ``self`` is a subclass of ``other``).
        """
        if isinstance(self, type(other)):
            first, second = (self, other) if not reverse else (other, self)
            return self.__class__(
                *[op(getattr(first, c), getattr(second, c)) for c in self.components]
            )
        else:
            try:
                self_cartesian = self.to_cartesian(other)
            except TypeError:
                return NotImplemented

            return other._combine_operation(op, self_cartesian, not reverse)

    def __sub__(self, other):
        # avoid "differential - representation".
        if isinstance(other, BaseRepresentation):
            return NotImplemented
        return super().__sub__(other)

    def norm(self, base=None):
        """Vector norm.

        The norm is the standard Frobenius norm, i.e., the square root of the
        sum of the squares of all components with non-angular units.

        Parameters
        ----------
        base : instance of ``self.base_representation``
            Base relative to which the differentials are defined. This is
            required to calculate the physical size of the differential for
            all but Cartesian differentials or radial differentials.

        Returns
        -------
        norm : `astropy.units.Quantity`
            Vector norm, with the same shape as the representation.
        """
        from .cartesian import CartesianDifferential

        # RadialDifferential overrides this function, so there is no handling here
        if not isinstance(self, CartesianDifferential) and base is None:
            raise ValueError(
                "`base` must be provided to calculate the norm of a"
                f" {type(self).__name__}"
            )
        return self.to_cartesian(base).norm()
