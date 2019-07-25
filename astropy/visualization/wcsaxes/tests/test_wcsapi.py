# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
import numpy as np

import matplotlib.pyplot as plt

from astropy.wcs.wcsapi import BaseLowLevelWCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity
from astropy.tests.image_tests import IMAGE_REFERENCE_DIR


class LowLevelWCS5D(BaseLowLevelWCS):

    @property
    def pixel_n_dim(self):
        return 2

    @property
    def world_n_dim(self):
        return 5

    @property
    def world_axis_physical_types(self):
        return ['em.freq', 'time', 'pos.eq.ra', 'pos.eq.dec', 'phys.polarization.stokes']

    @property
    def world_axis_units(self):
        return ['Hz', 'day', 'deg', 'deg', '']

    def pixel_to_world_values(self, *pixel_arrays):
        pixel_arrays = (list(pixel_arrays) * 3)[:-1]  # make list have 5 elements
        return [np.asarray(pix) * scale for pix, scale in zip(pixel_arrays, [10, 0.2, 0.4, 0.39, 2])]

    def array_index_to_world_values(self, *index_arrays):
        return self.pixel_to_world_values(index_arrays[::-1])[::-1]

    def world_to_pixel_values(self, *world_arrays):
        world_arrays = world_arrays[:2]  # make list have 2 elements
        return [np.asarray(world) / scale for world, scale in zip(world_arrays, [10, 0.2])]

    def world_to_array_index_values(self, *world_arrays):
        return np.round(self.world_to_array_index_values(world_arrays[::-1])[::-1]).astype(int)

    @property
    def world_axis_object_components(self):
        return [('freq', 0, 'value'),
                ('time', 0, 'mjd'),
                ('celestial', 0, 'spherical.lon.degree')
                ('celestial', 1, 'spherical.lat.degree'),
                ('stokes', 0, 'value')]

    @property
    def world_axis_object_classes(self):
        return {'celestial': (SkyCoord, (), {'unit': 'deg'}),
                'time': (Time, (), {'format': 'mjd'}),
                'freq': (Quantity, (), {'unit': 'Hz'}),
                'stokes': (Quantity, (), {'unit': 'one'})}


class TestWCSAPI:

    @pytest.mark.remote_data(source='astropy')
    @pytest.mark.mpl_image_compare(baseline_dir=IMAGE_REFERENCE_DIR,
                                   tolerance=0, style={})
    def test_wcsapi_5d(self):
        # Test for plotting image and also setting values of ticks
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=LowLevelWCS5D())
        ax.set_xlim(-0.5, 148.5)
        ax.set_ylim(-0.5, 148.5)
        return fig
