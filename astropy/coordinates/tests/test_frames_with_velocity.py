# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ... import units as u
from ..builtin_frames import ICRS, Galactic, Galactocentric

def test_api():
    # transform observed Barycentric velocities to full-space Galactocentric
    gc_frame = Galactocentric()
    icrs = ICRS(ra=151.*u.deg, dec=-16*u.deg, distance=101*u.pc,
                pm_ra_cosdec=21*u.mas/u.yr, pm_dec=-71*u.mas/u.yr,
                radial_velocity=71*u.km/u.s)
    icrs.transform_to(gc_frame)

    # transform a set of ICRS proper motions to Galactic
    icrs = ICRS(ra=151.*u.deg, dec=-16*u.deg,
                pm_ra_cosdec=21*u.mas/u.yr, pm_dec=-71*u.mas/u.yr)
    icrs.transform_to(Galactic)

    # transform a Barycentric RV to a GSR RV
    icrs = ICRS(ra=151.*u.deg, dec=-16*u.deg, distance=1.*u.pc,
                pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr,
                radial_velocity=71*u.km/u.s)
    icrs.transform_to(Galactocentric)

def test_all_arg_options():
    # I think this is a list of all possible valid combinations of arguments.
    # Here we do a simple thing and just verify that passing them in, we have
    # access to the relevant attributes from the resulting object
    all_kwargs = []

    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg)]
    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc)]

    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg,
                        pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr)]
    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
                        pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr)]

    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg,
                        radial_velocity=105.7*u.km/u.s)]
    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
                        radial_velocity=105.7*u.km/u.s)]

    all_kwargs += [dict(ra=37.4*u.deg, dec=-55.8*u.deg, distance=150*u.pc,
                        pm_ra_cosdec=-21.2*u.mas/u.yr, pm_dec=17.1*u.mas/u.yr,
                        radial_velocity=105.7*u.km/u.s)]

    for i,kwargs in enumerate(all_kwargs):
        icrs = ICRS(**kwargs)

        for k in kwargs:
            getattr(icrs, k)
