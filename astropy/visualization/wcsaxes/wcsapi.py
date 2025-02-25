# Functions/classes for WCSAxes related to APE14 WCSes
from __future__ import annotations

from contextlib import contextmanager

import numpy as np

from astropy import units as u
from astropy.coordinates import ICRS, BaseCoordinateFrame, SkyCoord
from astropy.wcs import WCS
from astropy.wcs.wcsapi import SlicedLowLevelWCS

from .frame import EllipticalFrame, RectangularFrame, RectangularFrame1D
from .transforms import CurvedTransform

__all__ = [
    "WCSPixel2WorldTransform",
    "WCSWorld2PixelTransform",
    "custom_ucd_coord_meta_mapping",
    "transform_coord_meta_from_wcs",
]

IDENTITY = WCS(naxis=2)
IDENTITY.wcs.ctype = ["X", "Y"]
IDENTITY.wcs.crval = [0.0, 0.0]
IDENTITY.wcs.crpix = [1.0, 1.0]
IDENTITY.wcs.cdelt = [1.0, 1.0]

UCD_COORD_META_MAPPING = {
    "lon": {"coord_type": "longitude"},
    "lat": {"coord_type": "latitude"},
    "ra": {"coord_type": "longitude", "format_unit": u.hourangle},
    "dec": {"coord_type": "latitude"},
    "alt": {"coord_type": "longitude"},
    "az": {"coord_type": "latitude"},
    "long": {"coord_type": "longitude"},
}

CUSTOM_UCD_COORD_META_MAPPING = {
    "pos.helioprojective.lon": {
        "coord_wrap": 180.0 * u.deg,
        "format_unit": u.arcsec,
        "coord_type": "longitude",
    },
    "pos.helioprojective.lat": {"format_unit": u.arcsec, "coord_type": "latitude"},
    "pos.heliographic.stonyhurst.lon": {
        "coord_wrap": 180.0 * u.deg,
        "format_unit": u.deg,
        "coord_type": "longitude",
    },
    "pos.heliographic.stonyhurst.lat": {"format_unit": u.deg, "coord_type": "latitude"},
    "pos.heliographic.carrington.lon": {
        "coord_wrap": 360.0 * u.deg,
        "format_unit": u.deg,
        "coord_type": "longitude",
    },
    "pos.heliographic.carrington.lat": {"format_unit": u.deg, "coord_type": "latitude"},
}


@contextmanager
def custom_ucd_coord_meta_mapping(mapping, *, overwrite=False):
    """
    A context manager that makes it possible to temporarily add new UCD+ to WCS coordinate
    plot metadata mappings.

    Parameters
    ----------
    mapping : dict
        A dictionary mapping a UCD to coordinate plot metadata.
        Note that custom UCD names have their "custom:" prefix stripped.
    overwrite : bool
        If `True` overwrite existing entries with ``mapping``.

    Examples
    --------
    >>> from matplotlib import pyplot as plt
    >>> from astropy.visualization.wcsaxes.wcsapi import custom_ucd_coord_meta_mapping
    >>> from astropy.wcs.wcsapi.fitswcs import custom_ctype_to_ucd_mapping
    >>> wcs = WCS(naxis=1)
    >>> wcs.wcs.ctype = ["eggs"]
    >>> wcs.wcs.cunit = ["deg"]
    >>> custom_mapping = {"eggs": "custom:pos.eggs"}
    >>> with custom_ctype_to_ucd_mapping(custom_mapping):
    ...     custom_meta = {
    ...         "pos.eggs": {
    ...             "coord_wrap": 360.0 * u.deg,
    ...             "format_unit": u.arcsec,
    ...             "coord_type": "longitude",
    ...         }
    ...     }
    ...     with custom_ucd_coord_meta_mapping(custom_meta):
    ...        fig = plt.figure()
    ...        ax = fig.add_subplot(111, projection=wcs)
    ...        ax.coords
    <CoordinatesMap with 1 world coordinates:
    <BLANKLINE>
      index       aliases           type   unit  wrap format_unit visible
                                                 deg
      ----- -------------------- --------- ---- ----- ----------- -------
          0 custom:pos.eggs eggs longitude  deg 360.0      arcsec     yes
    <BLANKLINE>
    >
    """
    added_keys = []
    overwritten = {}
    for k, v in mapping.items():
        k = k.removeprefix("custom:")
        if k in CUSTOM_UCD_COORD_META_MAPPING:
            if not overwrite:
                raise ValueError(f"UCD metadata mapping {k} already exists.")
            overwritten[k] = CUSTOM_UCD_COORD_META_MAPPING[k]
        else:
            added_keys.append(k)
        CUSTOM_UCD_COORD_META_MAPPING[k] = v

    try:
        yield
    finally:
        for k in added_keys:
            del CUSTOM_UCD_COORD_META_MAPPING[k]
        CUSTOM_UCD_COORD_META_MAPPING.update(overwritten)


def transform_coord_meta_from_wcs(wcs, frame_class, slices=None):
    if slices is not None:
        slices = tuple(slices)

    if wcs.pixel_n_dim > 2:
        if slices is None:
            raise ValueError(
                "WCS has more than 2 pixel dimensions, so 'slices' should be set"
            )
        elif len(slices) != wcs.pixel_n_dim:
            raise ValueError(
                "'slices' should have as many elements as WCS "
                f"has pixel dimensions (should be {wcs.pixel_n_dim})"
            )

    is_fits_wcs = isinstance(wcs, WCS) or (
        isinstance(wcs, SlicedLowLevelWCS) and isinstance(wcs._wcs, WCS)
    )

    coord_meta = {}
    coord_meta["name"] = []
    coord_meta["type"] = []
    coord_meta["wrap"] = []
    coord_meta["unit"] = []
    coord_meta["visible"] = []
    coord_meta["format_unit"] = []

    for idx in range(wcs.world_n_dim):
        axis_type = wcs.world_axis_physical_types[idx]
        axis_unit = u.Unit(wcs.world_axis_units[idx])
        coord_wrap = None
        format_unit = axis_unit
        coord_type = "scalar"

        dim_meta = {
            "coord_type": coord_type,
            "coord_wrap": coord_wrap,
            "format_unit": format_unit,
            "axis_unit": axis_unit,
        }

        if axis_type is not None:
            axis_type_split = axis_type.split(".")

            if len(axis_type_split):
                axis_type_split[0] = axis_type_split[0].replace("custom:", "")

            for ucd, meta in CUSTOM_UCD_COORD_META_MAPPING.items():
                if ucd in axis_type:
                    dim_meta.update(meta)
                    break
            else:
                for ucd, meta in UCD_COORD_META_MAPPING.items():
                    if ucd == axis_type_split[-1]:
                        dim_meta.update(meta)
                        break

        coord_meta["type"].append(dim_meta["coord_type"])
        coord_meta["wrap"].append(dim_meta["coord_wrap"])
        coord_meta["format_unit"].append(dim_meta["format_unit"])
        coord_meta["unit"].append(dim_meta["axis_unit"])

        # For FITS-WCS, for backward-compatibility, we need to make sure that we
        # provide aliases based on CTYPE for the name.
        if is_fits_wcs:
            name = []
            if isinstance(wcs, WCS):
                name.append(wcs.wcs.ctype[idx].lower())
                name.append(wcs.wcs.ctype[idx][:4].replace("-", "").lower())
            elif isinstance(wcs, SlicedLowLevelWCS):
                name.append(wcs._wcs.wcs.ctype[wcs._world_keep[idx]].lower())
                name.append(
                    wcs._wcs.wcs.ctype[wcs._world_keep[idx]][:4]
                    .replace("-", "")
                    .lower()
                )
            if name[0] == name[1]:
                name = name[0:1]
            if axis_type:
                if axis_type not in name:
                    name.insert(0, axis_type)
            if wcs.world_axis_names and wcs.world_axis_names[idx]:
                if wcs.world_axis_names[idx] not in name:
                    name.append(wcs.world_axis_names[idx])
            name = tuple(name) if len(name) > 1 else name[0]
        else:
            name = axis_type or ""
            if wcs.world_axis_names:
                name = (
                    (name, wcs.world_axis_names[idx])
                    if wcs.world_axis_names[idx]
                    else name
                )

        coord_meta["name"].append(name)

    coord_meta["default_axislabel_position"] = [""] * wcs.world_n_dim
    coord_meta["default_ticklabel_position"] = [""] * wcs.world_n_dim
    coord_meta["default_ticks_position"] = [""] * wcs.world_n_dim
    # If the world axis has a name use it, else display the world axis physical type.
    fallback_labels = [
        name[0] if isinstance(name, (list, tuple)) else name
        for name in coord_meta["name"]
    ]
    coord_meta["default_axis_label"] = [
        wcs.world_axis_names[i] or fallback_label
        for i, fallback_label in enumerate(fallback_labels)
    ]

    transform_wcs, invert_xy, world_map = apply_slices(wcs, slices)

    transform = WCSPixel2WorldTransform(transform_wcs, invert_xy=invert_xy)

    for i in range(len(coord_meta["type"])):
        coord_meta["visible"].append(i in world_map)

    inv_all_corr = [False] * wcs.world_n_dim
    m = transform_wcs.axis_correlation_matrix.copy()
    if invert_xy:
        inv_all_corr = np.all(m, axis=1)
        m = m[:, ::-1]

    if frame_class in (RectangularFrame, RectangularFrame1D):
        for index in world_map:
            coord_meta["default_axislabel_position"][index] = "#"
            coord_meta["default_ticklabel_position"][index] = "#"
            coord_meta["default_ticks_position"][index] = "#"

        # In the special and common case where the frame is rectangular and we
        # are dealing with a 2-d WCS (after slicing) for RectangularFrame or a
        # 1-d WCS for RectangularFrame1D, we show all ticks on all axes.
        if (frame_class is RectangularFrame and len(world_map) == 2) or (
            frame_class is RectangularFrame1D and len(world_map) == 1
        ):
            for index in world_map:
                coord_meta["default_ticks_position"][index] = frame_class.spine_names

    elif frame_class is EllipticalFrame:
        if "longitude" in coord_meta["type"]:
            lon_idx = coord_meta["type"].index("longitude")
            coord_meta["default_axislabel_position"][lon_idx] = "h"
            coord_meta["default_ticklabel_position"][lon_idx] = "h"
            coord_meta["default_ticks_position"][lon_idx] = "h"

        if "latitude" in coord_meta["type"]:
            lat_idx = coord_meta["type"].index("latitude")
            coord_meta["default_axislabel_position"][lat_idx] = "c"
            coord_meta["default_ticklabel_position"][lat_idx] = "c"
            coord_meta["default_ticks_position"][lat_idx] = "c"

    else:
        for index in range(len(coord_meta["type"])):
            if index in world_map:
                coord_meta["default_axislabel_position"][index] = (
                    frame_class.spine_names
                )
                coord_meta["default_ticklabel_position"][index] = (
                    frame_class.spine_names
                )
                coord_meta["default_ticks_position"][index] = frame_class.spine_names

    return transform, coord_meta


def apply_slices(wcs, slices):
    """
    Take the input WCS and slices and return a sliced WCS for the transform and
    a mapping of world axes in the sliced WCS to the input WCS.
    """
    if isinstance(wcs, SlicedLowLevelWCS):
        world_keep = list(wcs._world_keep)
    else:
        world_keep = list(range(wcs.world_n_dim))

    # world_map is the index of the world axis in the input WCS for a given
    # axis in the transform_wcs
    world_map = list(range(wcs.world_n_dim))
    transform_wcs = wcs
    invert_xy = False
    if slices is not None:
        wcs_slice = list(slices)
        wcs_slice[wcs_slice.index("x")] = slice(None)
        if "y" in slices:
            wcs_slice[wcs_slice.index("y")] = slice(None)
            invert_xy = slices.index("x") > slices.index("y")

        transform_wcs = SlicedLowLevelWCS(wcs, wcs_slice[::-1])
        world_map = tuple(world_keep.index(i) for i in transform_wcs._world_keep)

    return transform_wcs, invert_xy, world_map


def wcsapi_to_celestial_frame(wcs):
    for cls, _, kwargs, *_ in wcs.world_axis_object_classes.values():
        if issubclass(cls, SkyCoord):
            return kwargs.get("frame", ICRS())
        elif issubclass(cls, BaseCoordinateFrame):
            return cls(**kwargs)


class WCSWorld2PixelTransform(CurvedTransform):
    """
    WCS transformation from world to pixel coordinates.
    """

    has_inverse = True
    frame_in = None

    def __init__(self, wcs, invert_xy=False):
        super().__init__()

        if wcs.pixel_n_dim > 2:
            raise ValueError("Only pixel_n_dim =< 2 is supported")

        self.wcs = wcs
        self.invert_xy = invert_xy

        self.frame_in = wcsapi_to_celestial_frame(wcs)

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and self.wcs is other.wcs
            and self.invert_xy == other.invert_xy
        )

    @property
    def input_dims(self):
        return self.wcs.world_n_dim

    def transform(self, world):
        # Convert to a list of arrays
        world = list(world.T)

        if len(world) != self.wcs.world_n_dim:
            raise ValueError(
                f"Expected {self.wcs.world_n_dim} world coordinates, got {len(world)} "
            )

        if len(world[0]) == 0:
            pixel = np.zeros((0, 2))
        else:
            pixel = self.wcs.world_to_pixel_values(*world)

        if self.invert_xy:
            pixel = pixel[::-1]

        return np.array(pixel).T

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform.
        """
        return WCSPixel2WorldTransform(self.wcs, invert_xy=self.invert_xy)


class WCSPixel2WorldTransform(CurvedTransform):
    """
    WCS transformation from pixel to world coordinates.
    """

    has_inverse = True

    def __init__(self, wcs, invert_xy=False):
        super().__init__()

        if wcs.pixel_n_dim > 2:
            raise ValueError("Only pixel_n_dim =< 2 is supported")

        self.wcs = wcs
        self.invert_xy = invert_xy

        self.frame_out = wcsapi_to_celestial_frame(wcs)

    def __eq__(self, other):
        return (
            isinstance(other, type(self))
            and self.wcs is other.wcs
            and self.invert_xy == other.invert_xy
        )

    @property
    def output_dims(self):
        return self.wcs.world_n_dim

    def transform(self, pixel):
        # Convert to a list of arrays
        pixel = list(pixel.T)

        if len(pixel) != self.wcs.pixel_n_dim:
            raise ValueError(
                f"Expected {self.wcs.pixel_n_dim} world coordinates, got {len(pixel)} "
            )

        if self.invert_xy:
            pixel = pixel[::-1]

        if len(pixel[0]) == 0:
            world = np.zeros((0, self.wcs.world_n_dim))
        else:
            world = self.wcs.pixel_to_world_values(*pixel)

        if self.wcs.world_n_dim == 1:
            world = [world]

        return np.array(world).T

    transform_non_affine = transform

    def inverted(self):
        """
        Return the inverse of the transform.
        """
        return WCSWorld2PixelTransform(self.wcs, invert_xy=self.invert_xy)
