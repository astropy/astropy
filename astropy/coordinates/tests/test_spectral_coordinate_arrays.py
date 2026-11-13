# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Array metadata regressions for astropy/astropy#19560."""

import numpy as np
import pytest

import astropy.units as u
from astropy.constants import c
from astropy.coordinates import ICRS, CartesianDifferential, CartesianRepresentation
from astropy.coordinates.spectral_coordinate import NoVelocityWarning, SpectralCoord
from astropy.tests.helper import assert_quantity_allclose


def coordinates(shape=(2,)):
    n = np.prod(shape)
    xyz = np.arange(1, 3 * n + 1).reshape((*shape, 3)) * u.Mpc
    target = ICRS(
        CartesianRepresentation(
            xyz,
            xyz_axis=-1,
            differentials=CartesianDifferential(
                xyz.value * 10 * u.km / u.s, xyz_axis=-1
            ),
        )
    )
    observer = ICRS(
        CartesianRepresentation(
            [0, 0, 0] * u.kpc,
            differentials=CartesianDifferential([0, 0, 0] * u.km / u.s),
        )
    )
    return SpectralCoord(
        target.radial_velocity,
        doppler_convention="radio",
        doppler_rest=1.4 * u.GHz,
        observer=observer,
        target=target,
    )


@pytest.mark.parametrize(
    "operation",
    [
        lambda a: a[0],
        lambda a: a[:, 1:],
        lambda a: a[[1, 0]],
        lambda a: (
            a[a.value > a.value.mean()] if hasattr(a, "value") else a[a > a.mean()]
        ),
        lambda a: a.reshape(3, 2),
        lambda a: a.ravel(),
        lambda a: a.flatten(),
        lambda a: a.T,
        lambda a: a.swapaxes(0, 1),
        lambda a: a.diagonal(),
        lambda a: a.take([0, 2], axis=1),
        lambda a: a.repeat(2, axis=0),
        lambda a: a.compress([True, False], axis=0),
        lambda a: np.flip(a, axis=1),
        lambda a: np.roll(a, 1, axis=0),
        lambda a: np.expand_dims(a, axis=1),
        lambda a: np.moveaxis(a, 0, -1),
    ],
)
def test_metadata_shape_operations(operation):
    sc = coordinates((2, 3))
    expected = operation(sc.radial_velocity)
    result = operation(sc)
    assert result.shape == expected.shape
    assert result.target.shape == result.shape
    assert result.observer.isscalar
    assert_quantity_allclose(result.radial_velocity, expected)
    assert_quantity_allclose(result.quantity, operation(sc.quantity))


def test_broadcast_metadata_and_shape_assignment():
    sc = SpectralCoord(
        np.arange(6).reshape(2, 3) * u.GHz, radial_velocity=np.arange(3) * u.km / u.s
    )
    assert_quantity_allclose(sc[1].radial_velocity, sc.radial_velocity)
    expected = np.broadcast_to(sc.radial_velocity, sc.shape, subok=True).T
    result = sc.T
    assert_quantity_allclose(result.radial_velocity, expected)
    sc.shape = (3, 2)
    assert_quantity_allclose(
        sc.radial_velocity, np.array([0, 1, 2, 0, 1, 2]).reshape(3, 2) * u.km / u.s
    )
    assert sc.shape == sc.radial_velocity.shape


@pytest.mark.parametrize(
    "function, axis",
    [
        (np.concatenate, 0),
        (np.concatenate, 1),
        (np.concatenate, -1),
        (np.concatenate, None),
        (np.stack, 0),
        (np.stack, 1),
        (np.stack, -1),
    ],
)
def test_combine_particle_coordinates(function, axis):
    sc = coordinates((2, 3))
    result = function([sc, sc], axis=axis)
    assert_quantity_allclose(
        result.quantity, function([sc.quantity, sc.quantity], axis=axis)
    )
    assert_quantity_allclose(
        result.radial_velocity, function([sc.radial_velocity] * 2, axis=axis)
    )
    assert result.target.shape == result.shape
    assert result.observer is sc.observer
    assert result.doppler_rest is sc.doppler_rest


def test_combine_scalar_metadata_spectral_segments():
    sc = coordinates()[0].replicate(value=np.array([1, 2]) * u.GHz)
    other = sc.replicate(value=np.array([3000, 4000]) * u.MHz)
    result = np.concatenate([sc, other], dtype=np.float32)
    assert result.dtype == np.float32
    assert_quantity_allclose(result.quantity, [1, 2, 3, 4] * u.GHz)
    assert result.target.isscalar
    assert result.observer.isscalar
    assert result.radial_velocity.isscalar


def test_combine_distinct_scalar_radial_velocities():
    a = SpectralCoord([1, 2] * u.GHz, radial_velocity=10 * u.km / u.s)
    b = SpectralCoord([3, 4] * u.GHz, radial_velocity=20 * u.km / u.s)
    result = np.concatenate([a, b])
    assert_quantity_allclose(result.radial_velocity, [10, 10, 20, 20] * u.km / u.s)
    assert_quantity_allclose(result[2:].radial_velocity, [20, 20] * u.km / u.s)


def test_copy_metadata():
    sc = coordinates()
    result = sc.copy()
    result.target.data.x[0] += 1 * u.Mpc
    assert result.target.data.x[0] != sc.target.data.x[0]
    assert result.observer is not sc.observer
    assert not np.shares_memory(result, sc)


@pytest.mark.parametrize(
    "target_shift, observer_shift",
    [
        (10 * u.km / u.s, None),
        (None, 10 * u.km / u.s),
        (np.array([10, 20]) * u.km / u.s, None),
        (np.array([10, 20]) * u.km / u.s, np.array([3, 5]) * u.km / u.s),
    ],
)
def test_array_velocity_shift(target_shift, observer_shift):
    sc = coordinates()
    before = sc.radial_velocity.copy()
    result = sc.with_radial_velocity_shift(target_shift, observer_shift)
    dv = (0 * u.km / u.s if target_shift is None else target_shift) - (
        0 * u.km / u.s if observer_shift is None else observer_shift
    )
    assert_quantity_allclose(result.radial_velocity, before + dv)
    beta = (dv / c).to_value(u.one)
    factor = np.sqrt((1 + beta) / (1 - beta))
    assert_quantity_allclose(result.to(u.Hz).quantity, sc.to(u.Hz).quantity / factor)
    assert_quantity_allclose(sc.radial_velocity, before)


def test_array_to_rest():
    sc = coordinates()
    result = sc.to_rest()
    assert_quantity_allclose(
        result.radial_velocity, np.zeros(2) * u.km / u.s, atol=1e-10 * u.km / u.s
    )
    beta = (sc.radial_velocity / c).to_value(u.one)
    assert_quantity_allclose(
        result.to(u.Hz).quantity,
        sc.to(u.Hz).quantity * np.sqrt((1 + beta) / (1 - beta)),
    )


def test_array_without_velocity():
    target = ICRS(ra=[1, 2] * u.deg, dec=[3, 4] * u.deg, distance=[1, 2] * u.kpc)
    with pytest.warns(NoVelocityWarning):
        sc = SpectralCoord([1, 2] * u.GHz, target=target)
    assert sc.target.velocity.shape == (2,)
    assert_quantity_allclose(sc.target.velocity.d_xyz, np.zeros((3, 2)) * u.km / u.s)


def test_combine_incompatible_metadata():
    sc = coordinates()
    with pytest.raises(TypeError, match="All inputs"):
        np.concatenate([sc, sc.quantity])
    with pytest.raises(ValueError, match="Doppler convention"):
        np.concatenate([sc, sc.to(u.Hz, doppler_convention="optical")])
    other = SpectralCoord(
        [1, 2] * u.GHz, doppler_convention="radio", doppler_rest=1.4 * u.GHz
    )
    with pytest.raises(ValueError, match="observer"):
        np.concatenate([sc.to(u.GHz), other])
    with pytest.raises(NotImplementedError, match="out"):
        np.concatenate([sc, sc], out=np.empty(4) * sc.unit)


def test_shape_assignment_failure_is_atomic():
    sc = coordinates((2, 3)).T
    original_shape = sc.shape
    original_values = sc.quantity.copy()
    original_rv = sc.radial_velocity.copy()
    with pytest.raises(AttributeError, match="shape"):
        sc.shape = (6,)
    assert sc.shape == original_shape
    assert_quantity_allclose(sc.quantity, original_values)
    assert_quantity_allclose(sc.radial_velocity, original_rv)


def test_empty_selection_combination():
    sc = coordinates()
    empty = sc[:0]
    assert empty.target.shape == empty.shape == (0,)
    result = np.concatenate([empty, sc, empty])
    assert_quantity_allclose(result.quantity, sc.quantity)
    assert_quantity_allclose(result.radial_velocity, sc.radial_velocity)


@pytest.mark.parametrize("axis", [0, 1, -1, None])
def test_sort_metadata(axis):
    sc = coordinates((2, 3))[:, ::-1]
    result = np.sort(sc, axis=axis)
    expected = np.sort(sc.radial_velocity, axis=axis)
    assert_quantity_allclose(result.radial_velocity, expected)
    assert_quantity_allclose(result.quantity, np.sort(sc.quantity, axis=axis))
    assert sc.radial_velocity[0, 0] > sc.radial_velocity[0, 1]


def test_inplace_sort_and_partition_reject_array_metadata():
    sc = coordinates()[::-1]
    before = sc.quantity.copy()
    with pytest.raises(NotImplementedError, match="In-place sorting"):
        sc.sort()
    for operation in (lambda: sc.partition(0), lambda: np.partition(sc, 0)):
        with pytest.raises(NotImplementedError, match="Partitioning"):
            operation()
    assert_quantity_allclose(sc.quantity, before)


def test_slice_metadata_keeps_view_semantics():
    sc = SpectralCoord([1, 2] * u.GHz, radial_velocity=[10, 20] * u.km / u.s)
    sliced = sc[:1]
    sliced.radial_velocity[0] = 77 * u.km / u.s
    assert_quantity_allclose(sc.radial_velocity, [77, 20] * u.km / u.s)


@pytest.mark.parametrize("order", ["C", "F", "K"])
def test_numpy_copy_positional_order(order):
    sc = coordinates((2, 3))
    result = np.copy(sc, order)
    assert_quantity_allclose(result.radial_velocity, sc.radial_velocity)
    assert not np.shares_memory(result, sc)
    assert not np.shares_memory(result.target.data.x, sc.target.data.x)

