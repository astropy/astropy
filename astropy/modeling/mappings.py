"""
Special models useful for complex compound models where control is needed over
which outputs from a source model are mapped to which inputs of a target model.
"""
# pylint: disable=invalid-name

from .core import FittableModel


__all__ = ['Mapping', 'Identity']


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
