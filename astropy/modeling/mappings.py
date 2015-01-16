"""
Special models useful for complex compound models where control is needed over
which outputs from a source model are mapped to which inputs of a target model.
"""

from .core import Model
from .functional_models import Shift, Scale
from .polynomial import Polynomial2D


__all__ = ['Mapping', 'Identity']


class Mapping(Model):
    """
    Allows inputs to be reordered, duplicated or dropped.

    Parameters
    ----------
    mapping : tuple
        Integers representing indices of the inputs.
    n_inputs : int
        Number of inputs; if None(default) then ``max(mapping)`` + 1 is used.
    name : str
        Name for this model.
    meta : dict-like
        Metadata

    Raises
    ------
    TypeError
        Raised when number of inputs is less that ``max(mapping)``.

    Examples
    --------
    >>> poly1 = Polynomial2D(1, c0_0=1, c1_0=2, c0_1=3)
    >>> poly2 = Polynomial2D(1, c0_0=1, c1_0=2.4, c0_1=2.1)
    >>> model = (Shift(1) & Shift(2)) | Mapping((0, 1, 0, 1)) | (poly1 & poly2)
    >>> model(1, 2)
    (17.0, 14.2)

    """
    def __init__(self, mapping, n_inputs=None, name=None, meta=None):
        if n_inputs is None:
            self._inputs = tuple('x' + str(idx)
                                 for idx in range(max(mapping) + 1))
        else:
            self._inputs = tuple('x' + str(idx)
                                 for idx in range(n_inputs))
        self._outputs = tuple('x' + str(idx) for idx in range(len(mapping)))
        self._mapping = mapping
        super(Mapping, self).__init__(name=name, meta=meta)

    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs

    @property
    def mapping(self):
        return self._mapping

    def __repr__(self):
        if self.name is None:
            return '<Mapping({0})>'.format(self.mapping)
        else:
            return '<Mapping({0}, name={1})>'.format(self.mapping, self.name)

    def evaluate(self, *args):
        if len(args) != self.n_inputs:
            raise TypeError('{0} expects {1} inputs; got {2}'.format(
                self.name, self.n_inputs, len(args)))

        result = tuple(args[idx] for idx in self._mapping)

        if self.n_outputs == 1:
            return result[0]

        return result

    @property
    def inverse(self):
        try:
            mapping = tuple(self.mapping.index(idx)
                            for idx in range(self.n_inputs))
        except ValueError:
            raise NotImplementedError(
                "Mappings such as {0} that drop one or more of their inputs "
                "are not invertible at this time.".format(self.mapping))

        inv = self.__class__(mapping)
        inv._inputs = self._outputs
        inv._outputs = self._inputs
        return inv


class Identity(Mapping):
    """
    Returns inputs unchanged.

    This class is useful in compound models when some of the inputs must be
    passed unchanged to the next model.

    Parameters
    ----------
    n_inputs : int
        Specifies how many of the inputs will be returned.
    name : str
        User supplied model name.
    meta : dict-like
        Metadata

    Examples
    --------
    >>> # Transform (x, y) by a shift in x, followed by scaling the two inputs.
    >>> model = (Shift(1) & Identity(1)) | Scale(1.2) & Scale(2)
    >>> model(1,1)
    (2.4, 2.0)
    >>> model.inverse(2.4, 2)
    (1.0, 1.0)

    """
    def __init__(self, n_inputs, name=None, meta=None):
        mapping = tuple(range(n_inputs))
        super(Identity, self).__init__(mapping, name=name, meta=meta)

    def __repr__(self):
        if self.name is None:
            return '<Identity({0})>'.format(self.n_inputs)
        else:
            return '<Identity({0}, name={1})>'.format(self.n_inputs, self.name)

    @property
    def inverse(self):
        return self
