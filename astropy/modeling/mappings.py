"""
Special models useful for complex compound models where control is needed over
which outputs from a source model are mapped to which inputs of a target model.
"""
# pylint: disable=invalid-name

from .core import FittableModel, Model
from astropy.units import Quantity


__all__ = ['Mapping', 'Identity', 'UnitsMapping']


class Mapping(FittableModel):
    """
    Allows inputs to be reordered, duplicated or dropped.

    Parameters
    ----------
    mapping : tuple
        A tuple of integers representing indices of the inputs to this model
        to return and in what order to return them.  See
        :ref:`compound-model-mappings` for more details.
    n_inputs : int
        Number of inputs; if `None` (default) then ``max(mapping) + 1`` is
        used (i.e. the highest input index used in the mapping).
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).
    meta : dict-like
        Free-form metadata to associate with this model.

    Raises
    ------
    TypeError
        Raised when number of inputs is less that ``max(mapping)``.

    Examples
    --------

    >>> from astropy.modeling.models import Polynomial2D, Shift, Mapping
    >>> poly1 = Polynomial2D(1, c0_0=1, c1_0=2, c0_1=3)
    >>> poly2 = Polynomial2D(1, c0_0=1, c1_0=2.4, c0_1=2.1)
    >>> model = (Shift(1) & Shift(2)) | Mapping((0, 1, 0, 1)) | (poly1 & poly2)
    >>> model(1, 2)  # doctest: +FLOAT_CMP
    (17.0, 14.2)
    """
    linear = True  # FittableModel is non-linear by default

    def __init__(self, mapping, n_inputs=None, name=None, meta=None):
        self._inputs = ()
        self._outputs = ()
        if n_inputs is None:
            self._n_inputs = max(mapping) + 1
        else:
            self._n_inputs = n_inputs

        self._n_outputs = len(mapping)
        super().__init__(name=name, meta=meta)

        self.inputs = tuple('x' + str(idx) for idx in range(self._n_inputs))
        self.outputs = tuple('x' + str(idx) for idx in range(self._n_outputs))

        self._mapping = mapping
        self._input_units_strict = {key: False for key in self._inputs}
        self._input_units_allow_dimensionless = {key: False for key in self._inputs}

    @property
    def n_inputs(self):
        return self._n_inputs

    @property
    def n_outputs(self):
        return self._n_outputs

    @property
    def mapping(self):
        """Integers representing indices of the inputs."""
        return self._mapping

    def __repr__(self):
        if self.name is None:
            return f'<Mapping({self.mapping})>'
        return f'<Mapping({self.mapping}, name={self.name!r})>'

    def evaluate(self, *args):
        if len(args) != self.n_inputs:
            name = self.name if self.name is not None else "Mapping"

            raise TypeError('{} expects {} inputs; got {}'.format(
                name, self.n_inputs, len(args)))

        result = tuple(args[idx] for idx in self._mapping)

        if self.n_outputs == 1:
            return result[0]

        return result

    @property
    def inverse(self):
        """
        A `Mapping` representing the inverse of the current mapping.

        Raises
        ------
        `NotImplementedError`
            An inverse does no exist on mappings that drop some of its inputs
            (there is then no way to reconstruct the inputs that were dropped).
        """

        try:
            mapping = tuple(self.mapping.index(idx)
                            for idx in range(self.n_inputs))
        except ValueError:
            raise NotImplementedError(
                "Mappings such as {} that drop one or more of their inputs "
                "are not invertible at this time.".format(self.mapping))

        inv = self.__class__(mapping)
        inv._inputs = self._outputs
        inv._outputs = self._inputs
        inv._n_inputs = len(inv._inputs)
        inv._n_outputs = len(inv._outputs)
        return inv


class Identity(Mapping):
    """
    Returns inputs unchanged.

    This class is useful in compound models when some of the inputs must be
    passed unchanged to the next model.

    Parameters
    ----------
    n_inputs : int
        Specifies the number of inputs this identity model accepts.
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).
    meta : dict-like
        Free-form metadata to associate with this model.

    Examples
    --------

    Transform ``(x, y)`` by a shift in x, followed by scaling the two inputs::

        >>> from astropy.modeling.models import (Polynomial1D, Shift, Scale,
        ...                                      Identity)
        >>> model = (Shift(1) & Identity(1)) | Scale(1.2) & Scale(2)
        >>> model(1,1)  # doctest: +FLOAT_CMP
        (2.4, 2.0)
        >>> model.inverse(2.4, 2) # doctest: +FLOAT_CMP
        (1.0, 1.0)
    """
    linear = True  # FittableModel is non-linear by default

    def __init__(self, n_inputs, name=None, meta=None):
        mapping = tuple(range(n_inputs))
        super().__init__(mapping, name=name, meta=meta)

    def __repr__(self):
        if self.name is None:
            return f'<Identity({self.n_inputs})>'
        return f'<Identity({self.n_inputs}, name={self.name!r})>'

    @property
    def inverse(self):
        """
        The inverse transformation.

        In this case of `Identity`, ``self.inverse is self``.
        """

        return self


class UnitsMapping(Model):
    """
    Mapper that operates on the units of the input, first converting to
    canonical units, then assigning new units without further conversion.
    Used by Model.coerce_units to support units on otherwise unitless models
    such as Polynomial1D.

    Parameters
    ----------
    mapping : tuple
        A tuple of (input_unit, output_unit) pairs, one per input, matched to the
        inputs by position.  The first element of the each pair is the unit that
        the model will accept (specify ``dimensionless_unscaled``
        to accept dimensionless input).  The second element is the unit that the
        model will return.  Specify ``dimensionless_unscaled``
        to return dimensionless Quantity, and `None` to return raw values without
        Quantity.
    input_units_equivalencies : dict, optional
        Default equivalencies to apply to input values.  If set, this should be a
        dictionary where each key is a string that corresponds to one of the
        model inputs.
    input_units_allow_dimensionless : dict or bool, optional
        Allow dimensionless input. If this is True, input values to evaluate will
        gain the units specified in input_units. If this is a dictionary then it
        should map input name to a bool to allow dimensionless numbers for that
        input.
    name : str, optional
        A human-friendly name associated with this model instance
        (particularly useful for identifying the individual components of a
        compound model).
    meta : dict-like, optional
        Free-form metadata to associate with this model.

    Examples
    --------

    Wrapping a unitless model to require and convert units:

    >>> from astropy.modeling.models import Polynomial1D, UnitsMapping
    >>> from astropy import units as u
    >>> poly = Polynomial1D(1, c0=1, c1=2)
    >>> model = UnitsMapping(((u.m, None),)) | poly
    >>> model = model | UnitsMapping(((None, u.s),))
    >>> model(u.Quantity(10, u.m))  # doctest: +FLOAT_CMP
    <Quantity 21. s>
    >>> model(u.Quantity(1000, u.cm)) # doctest: +FLOAT_CMP
    <Quantity 21. s>
    >>> model(u.Quantity(10, u.cm)) # doctest: +FLOAT_CMP
    <Quantity 1.2 s>

    Wrapping a unitless model but still permitting unitless input:

    >>> from astropy.modeling.models import Polynomial1D, UnitsMapping
    >>> from astropy import units as u
    >>> poly = Polynomial1D(1, c0=1, c1=2)
    >>> model = UnitsMapping(((u.m, None),), input_units_allow_dimensionless=True) | poly
    >>> model = model | UnitsMapping(((None, u.s),))
    >>> model(u.Quantity(10, u.m))  # doctest: +FLOAT_CMP
    <Quantity 21. s>
    >>> model(10)  # doctest: +FLOAT_CMP
    <Quantity 21. s>
    """
    def __init__(
        self,
        mapping,
        input_units_equivalencies=None,
        input_units_allow_dimensionless=False,
        name=None,
        meta=None
    ):
        self._mapping = mapping

        none_mapping_count = len([m for m in mapping if m[-1] is None])
        if none_mapping_count > 0 and none_mapping_count != len(mapping):
            raise ValueError("If one return unit is None, then all must be None")

        # These attributes are read and handled by Model
        self._input_units_strict = True
        self.input_units_equivalencies = input_units_equivalencies
        self._input_units_allow_dimensionless = input_units_allow_dimensionless

        super().__init__(name=name, meta=meta)

        # Can't invoke this until after super().__init__, since
        # we need self.inputs and self.outputs to be populated.
        self._rebuild_units()

    def _rebuild_units(self):
        self._input_units = {input_name: input_unit for input_name, (input_unit, _) in zip(self.inputs, self.mapping)}

    @property
    def n_inputs(self):
        return len(self._mapping)

    @property
    def n_outputs(self):
        return len(self._mapping)

    @property
    def inputs(self):
        return super().inputs

    @inputs.setter
    def inputs(self, value):
        super(UnitsMapping, self.__class__).inputs.fset(self, value)
        self._rebuild_units()

    @property
    def outputs(self):
        return super().outputs

    @outputs.setter
    def outputs(self, value):
        super(UnitsMapping, self.__class__).outputs.fset(self, value)
        self._rebuild_units()

    @property
    def input_units(self):
        return self._input_units

    @property
    def mapping(self):
        return self._mapping

    def evaluate(self, *args):
        result = []
        for arg, (_, return_unit) in zip(args, self.mapping):
            if isinstance(arg, Quantity):
                value = arg.value
            else:
                value = arg
            if return_unit is None:
                result.append(value)
            else:
                result.append(Quantity(value, return_unit, subok=True))

        if self.n_outputs == 1:
            return result[0]
        else:
            return tuple(result)

    def __repr__(self):
        if self.name is None:
            return f"<UnitsMapping({self.mapping})>"
        else:
            return f"<UnitsMapping({self.mapping}, name={self.name!r})>"
