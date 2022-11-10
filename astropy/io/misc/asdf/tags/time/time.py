# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
from asdf.versioning import AsdfSpec
from numpy.testing import assert_array_equal

from astropy import time
from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.io.misc.asdf.types import AstropyAsdfType
from astropy.units import Quantity

__all__ = ["TimeType"]

_guessable_formats = {"iso", "byear", "jyear", "yday"}

_astropy_format_to_asdf_format = {
    "isot": "iso",
    "byear_str": "byear",
    "jyear_str": "jyear",
}


def _assert_earthlocation_equal(a, b):
    assert_array_equal(a.x, b.x)
    assert_array_equal(a.y, b.y)
    assert_array_equal(a.z, b.z)
    assert_array_equal(a.lat, b.lat)
    assert_array_equal(a.lon, b.lon)


class TimeType(AstropyAsdfType):
    name = "time/time"
    version = "1.1.0"
    supported_versions = ["1.0.0", AsdfSpec(">=1.1.0")]
    types = ["astropy.time.core.Time"]
    requires = ["astropy"]

    @classmethod
    def to_tree(cls, node, ctx):
        fmt = node.format

        if fmt == "byear":
            node = time.Time(node, format="byear_str")

        elif fmt == "jyear":
            node = time.Time(node, format="jyear_str")

        elif fmt in ("fits", "datetime", "plot_date"):
            node = time.Time(node, format="isot")

        fmt = node.format

        fmt = _astropy_format_to_asdf_format.get(fmt, fmt)

        guessable_format = fmt in _guessable_formats

        if node.scale == "utc" and guessable_format and node.isscalar:
            return node.value

        d = {"value": node.value}

        if not guessable_format:
            d["format"] = fmt

        if node.scale != "utc":
            d["scale"] = node.scale

        if node.location is not None:
            x, y, z = node.location.x, node.location.y, node.location.z
            # Preserve backwards compatibility for writing the old schema
            # This allows WCS to test backwards compatibility with old frames
            # This code does get tested in CI, but we don't run a coverage test
            if cls.version == "1.0.0":  # pragma: no cover
                unit = node.location.unit
                d["location"] = {"x": x.value, "y": y.value, "z": z.value, "unit": unit}
            else:
                d["location"] = {
                    # It seems like EarthLocations can be represented either in
                    # terms of Cartesian coordinates or latitude and longitude, so
                    # we rather arbitrarily choose the former for our representation
                    "x": x,
                    "y": y,
                    "z": z,
                }

        return d

    @classmethod
    def from_tree(cls, node, ctx):
        if isinstance(node, (str, list, np.ndarray)):
            t = time.Time(node)
            fmt = _astropy_format_to_asdf_format.get(t.format, t.format)
            if fmt not in _guessable_formats:
                raise ValueError(f"Invalid time '{node}'")
            return t

        value = node["value"]
        fmt = node.get("format")
        scale = node.get("scale")
        location = node.get("location")
        if location is not None:
            unit = location.get("unit", u.m)
            # This ensures that we can read the v.1.0.0 schema and convert it
            # to the new EarthLocation object, which expects Quantity components
            for comp in ["x", "y", "z"]:
                if not isinstance(location[comp], Quantity):
                    location[comp] = Quantity(location[comp], unit=unit)
            location = EarthLocation.from_geocentric(
                location["x"], location["y"], location["z"]
            )

        return time.Time(value, format=fmt, scale=scale, location=location)

    @classmethod
    def assert_equal(cls, old, new):
        assert old.format == new.format
        assert old.scale == new.scale
        if isinstance(old.location, EarthLocation):
            assert isinstance(new.location, EarthLocation)
            _assert_earthlocation_equal(old.location, new.location)
        else:
            assert old.location == new.location

        assert_array_equal(old, new)
