"""
Special models useful for complex compound models where control is needed over
which outputs from a source model are mapped to which inputs of a target model.
"""

from .core import Model


__all__ = ['Mapping', 'Identity']


class Mapping(Model):
    def __init__(self, mapping, n_inputs=None):
        if n_inputs is None:
            self._inputs = tuple('x' + str(idx)
                                 for idx in range(max(mapping) + 1))
        else:
            self._inputs = tuple('x' + str(idx)
                                 for idx in range(n_inputs))
        self._outputs = tuple('x' + str(idx) for idx in range(len(mapping)))
        self._mapping = mapping
        super(Mapping, self).__init__()

    @property
    def name(self):
        return 'Mapping({0})'.format(self.mapping)

    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs

    @property
    def mapping(self):
        return self._mapping

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
    def __init__(self, n_inputs):
        mapping = tuple(range(n_inputs))
        super(Identity, self).__init__(mapping, n_inputs=n_inputs)

    @property
    def name(self):
        return 'Identity({0})'.format(self.n_inputs)

    @property
    def inverse(self):
        return self
