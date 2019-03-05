"""
This is a utility for working out the correlation between inputs and outputs of
a model.

This is subtly different from the separability matrix, as if two inputs go
through a single model they are coupled irrespective of the mathematical
separability of that model.
"""


import numpy as np

from astropy.modeling import separable
from astropy.modeling.core import Model, _CompoundModel, ModelDefinitionError
from astropy.modeling.mappings import Mapping


__all__ = ['coupling_matrix']


def _coord_matrix(model, pos, noutp):
    """
    Create an array representing inputs and outputs of a simple model.

    The array has a shape (noutp, model.n_inputs).

    Parameters
    ----------
    model : `astropy.modeling.Model`
        model
    pos : str
        Position of this model in the expression tree.
        One of ['left', 'right'].
    noutp : int
        Number of outputs of the compound model of which the input model
        is a left or right child.

    """
    if isinstance(model, Mapping):
        axes = []
        for i in model.mapping:
            axis = np.zeros((model.n_inputs,))
            axis[i] = 1
            axes.append(axis)
        m = np.vstack(axes)
        mat = np.zeros((noutp, model.n_inputs))
        if pos == 'left':
            mat[: model.n_outputs, :model.n_inputs] = m
        else:
            mat[-model.n_outputs:, -model.n_inputs:] = m
        return mat
    else:
        # this does not work for more than 2 coordinates
        mat = np.zeros((noutp, model.n_inputs))
        if pos == 'left':
            mat[:model.n_outputs, : model.n_inputs] = 1
        else:
            mat[-model.n_outputs:, -model.n_inputs:] = 1
    return mat


def _cdot(left, right):
    """
    Function corresponding to "|" operation.

    Parameters
    ----------
    left, right : `astropy.modeling.Model` or ndarray
        If input is of an array, it is the output of `coord_matrix`.

    Returns
    -------
    result : ndarray
        Result from this operation.
    """

    left, right = right, left

    def _n_inputs_outputs(input, position):
        """
        Return ``n_inputs``, ``n_outputs`` for a model or coord_matrix.
        """
        if isinstance(input, Model):
            coords = _coord_matrix(input, position, input.n_outputs)
        else:
            coords = input
        return coords

    cleft = _n_inputs_outputs(left, 'left')
    cright = _n_inputs_outputs(right, 'right')

    try:
        result = np.dot(cleft, cright)
    except ValueError:
        raise ModelDefinitionError(
            'Models cannot be combined with the "|" operator; '
            'left coord_matrix is {0}, right coord_matrix is {1}'.format(
                cright, cleft))
    return result


def _cstack(left, right):
    """
    Function corresponding to '&' operation.

    Parameters
    ----------
    left, right : `astropy.modeling.Model` or ndarray
        If input is of an array, it is the output of `coord_matrix`.

    Returns
    -------
    result : ndarray
        Result from this operation.

    """
    noutp = separable._compute_n_outputs(left, right)

    if isinstance(left, Model):
        cleft = _coord_matrix(left, 'left', noutp)
    else:
        cleft = np.zeros((noutp, left.shape[1]))
        cleft[: left.shape[0], : left.shape[1]] = left
    if isinstance(right, Model):
        cright = _coord_matrix(right, 'right', noutp)
    else:
        cright = np.zeros((noutp, right.shape[1]))
        cright[-right.shape[0]:, -right.shape[1]:] = 1

    return np.hstack([cleft, cright])


def _coupled(transform):
    """
    Calculate the separability of outputs.

    Parameters
    ----------
    transform : `astropy.modeling.Model`
        A transform (usually a compound model).

    Returns
    -------
    is_separable : ndarray of dtype np.bool
        An array of shape (transform.n_outputs,) of boolean type
        Each element represents the separablity of the corresponding output.
    """
    if isinstance(transform, _CompoundModel):
        is_separable = transform._tree.evaluate(_operators)
    elif isinstance(transform, Model):
        is_separable = _coord_matrix(transform, 'left', transform.n_outputs)
    return is_separable


def coupling_matrix(transform):
    """
    Compute the coupling between outputs and inputs.

    Parameters
    ----------
    transform : `~astropy.modeling.core.Model`
        A (compound) model.

    Returns
    -------
    coupling_matrix : ndarray
        A boolean correlation matrix of shape (n_outputs, n_inputs).
        Indicates the dependence of outputs on inputs. For completely
        independent outputs, the diagonal elements are True and
        off-diagonal elements are False.

    """
    if transform.n_inputs == 1 and transform.n_outputs > 1:
        return np.ones((transform.n_outputs, transform.n_inputs),
                       dtype=np.bool)
    separable_matrix = _coupled(transform)
    separable_matrix = np.where(separable_matrix != 0, True, False)
    return separable_matrix


_operators = separable._operators
# Maps modeling operators to a function computing and represents the
# relationship of axes as an array of 0-es and 1-s
_operators['|'] = _cdot
_operators['&'] = _cstack
