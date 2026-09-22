# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Tests for `~astropy.coordinates.get_frame_class`.

Covers:
- Resolution of registered frame name strings to frame classes.
- Acceptance of ``BaseCoordinateFrame`` subclasses (idempotent pass-through).
- Rejection of frame *instances* (must be a class, not an object).
- Rejection of unrecognised string names.
- Rejection of unrelated object types.
- Public import path ``astropy.coordinates.get_frame_class``.
- Backward-compatibility alias ``_get_frame_class``.
"""

import pytest

import astropy.coordinates as coord
from astropy.coordinates import (
    FK4,
    FK5,
    GCRS,
    ICRS,
    BaseCoordinateFrame,
    Galactic,
)
from astropy.coordinates.sky_coordinate_parsers import (
    _get_frame_class,
    get_frame_class,
)


# ---------------------------------------------------------------------------
# Helper: all frames expected to round-trip through their registered name
# ---------------------------------------------------------------------------
NAMED_FRAME_CASES = [
    ("icrs", ICRS),
    ("fk5", FK5),
    ("fk4", FK4),
    ("galactic", Galactic),
    ("gcrs", GCRS),
]


class TestGetFrameClassStringInput:
    """``get_frame_class`` called with a frame name string."""

    @pytest.mark.parametrize("name,expected_cls", NAMED_FRAME_CASES)
    def test_known_frame_names(self, name, expected_cls):
        """Registered name strings resolve to the correct class."""
        result = get_frame_class(name)
        assert result is expected_cls

    def test_unknown_frame_name_raises(self):
        """An unregistered string raises ``ValueError``."""
        with pytest.raises(ValueError, match="not a known coordinate frame"):
            get_frame_class("completely_unknown_frame_xyz")

    def test_empty_string_raises(self):
        """An empty string is not a registered name and must raise."""
        with pytest.raises(ValueError, match="not a known coordinate frame"):
            get_frame_class("")

    def test_case_sensitive(self):
        """Frame names are case-sensitive; 'ICRS' should not resolve."""
        with pytest.raises(ValueError, match="not a known coordinate frame"):
            get_frame_class("ICRS")


class TestGetFrameClassClassInput:
    """``get_frame_class`` called with a frame *class* (idempotent)."""

    @pytest.mark.parametrize("frame_cls", [ICRS, FK5, FK4, Galactic, GCRS])
    def test_frame_class_returns_itself(self, frame_cls):
        """A ``BaseCoordinateFrame`` subclass is returned unchanged."""
        result = get_frame_class(frame_cls)
        assert result is frame_cls

    def test_base_coordinate_frame_itself_is_accepted(self):
        """``BaseCoordinateFrame`` itself (the abstract base) is a valid class input."""
        result = get_frame_class(BaseCoordinateFrame)
        assert result is BaseCoordinateFrame

    def test_custom_frame_subclass(self):
        """A user-defined ``BaseCoordinateFrame`` subclass is accepted."""
        class MyCustomFrame(BaseCoordinateFrame):
            default_representation = coord.SphericalRepresentation

        result = get_frame_class(MyCustomFrame)
        assert result is MyCustomFrame


class TestGetFrameClassInvalidInput:
    """``get_frame_class`` must raise ``ValueError`` for invalid inputs."""

    def test_frame_instance_raises(self):
        """A frame *instance* (not class) must raise ``ValueError``."""
        instance = ICRS()
        with pytest.raises(ValueError, match="Coordinate frame must be a frame name or frame class"):
            get_frame_class(instance)

    def test_none_raises(self):
        """``None`` is not a valid frame specifier."""
        with pytest.raises(ValueError, match="Coordinate frame must be a frame name or frame class"):
            get_frame_class(None)

    def test_integer_raises(self):
        """An integer is not a valid frame specifier."""
        with pytest.raises(ValueError, match="Coordinate frame must be a frame name or frame class"):
            get_frame_class(42)

    def test_list_raises(self):
        """A list is not a valid frame specifier."""
        with pytest.raises(ValueError, match="Coordinate frame must be a frame name or frame class"):
            get_frame_class(["icrs"])

    def test_unrelated_class_raises(self):
        """A class unrelated to ``BaseCoordinateFrame`` must raise."""
        with pytest.raises(ValueError, match="Coordinate frame must be a frame name or frame class"):
            get_frame_class(int)


class TestGetFrameClassPublicAPI:
    """Verify the public import path works correctly."""

    def test_importable_from_astropy_coordinates(self):
        """``get_frame_class`` must be importable from ``astropy.coordinates``."""
        from astropy.coordinates import get_frame_class as _gfc  # noqa: F401
        assert callable(_gfc)

    def test_top_level_namespace_symbol(self):
        """``astropy.coordinates.get_frame_class`` is the correct public symbol."""
        assert hasattr(coord, "get_frame_class")
        assert coord.get_frame_class is get_frame_class

    def test_public_api_resolves_icrs(self):
        """End-to-end: the top-level public API resolves 'icrs' to ICRS."""
        result = coord.get_frame_class("icrs")
        assert result is ICRS

    def test_public_api_accepts_class(self):
        """End-to-end: the top-level public API accepts a frame class."""
        result = coord.get_frame_class(FK5)
        assert result is FK5


class TestGetFrameClassBackwardCompatAlias:
    """Verify the private backward-compatibility alias ``_get_frame_class`` still works."""

    def test_private_alias_is_same_function(self):
        """``_get_frame_class`` is the same object as ``get_frame_class``."""
        assert _get_frame_class is get_frame_class

    def test_private_alias_resolves_icrs(self):
        """``_get_frame_class('icrs')`` still works via the alias."""
        result = _get_frame_class("icrs")
        assert result is ICRS

    def test_private_alias_raises_on_invalid(self):
        """``_get_frame_class`` raises on invalid input, same as the public name."""
        with pytest.raises(ValueError):
            _get_frame_class("nonexistent_frame_xyz")


class TestGetFrameClassReturnType:
    """Verify the return value is always a class, never an instance."""

    @pytest.mark.parametrize("name,expected_cls", NAMED_FRAME_CASES)
    def test_returns_class_not_instance(self, name, expected_cls):
        """Return value must be a class (``type``), not an object."""
        result = get_frame_class(name)
        assert isinstance(result, type)
        assert issubclass(result, BaseCoordinateFrame)

    def test_class_can_be_instantiated(self):
        """The returned class is directly instantiable."""
        cls = get_frame_class("fk5")
        instance = cls()
        assert isinstance(instance, BaseCoordinateFrame)
