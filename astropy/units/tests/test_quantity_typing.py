# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test the Quantity class and related."""

from typing import Annotated, Any, Optional, Union, get_args, get_origin, get_type_hints

import numpy as np

from astropy import units as u


class TestQuantityTyping:
    """Test Quantity Typing Annotations."""

    def test_quantity_typing(self):
        """Test type hint creation from Quantity."""
        annot = u.Quantity[u.m]

        assert get_origin(annot) is Annotated
        assert get_args(annot) == (u.Quantity, u.m)

        # test usage
        def func(x: annot, y: str) -> u.Quantity[u.s]:
            return x, y

        annots = get_type_hints(func, include_extras=True)
        assert annots["x"] is annot
        assert annots["return"].__metadata__[0] == u.s

    def test_metadata_in_annotation(self):
        """Test Quantity annotation with added metadata."""
        multi_annot = u.Quantity[u.m, Any, np.dtype]

        def multi_func(x: multi_annot, y: str):
            return x, y

        annots = get_type_hints(multi_func, include_extras=True)
        assert annots["x"] == multi_annot

    def test_optional_and_annotated(self):
        """Test Quantity annotation in an Optional."""
        opt_annot = Optional[u.Quantity[u.m]]

        def opt_func(x: opt_annot, y: str):
            return x, y

        annots = get_type_hints(opt_func, include_extras=True)
        assert annots["x"] == opt_annot

    def test_union_and_annotated(self):
        """Test Quantity annotation in a Union."""
        # double Quantity[]
        union_annot1 = Union[u.Quantity[u.m], u.Quantity[u.s]]
        # one Quantity, one physical-type
        union_annot2 = Union[u.Quantity[u.m], u.Quantity["time"]]
        # one Quantity, one general type
        union_annot3 = Union[u.Quantity[u.m / u.s], float]

        def union_func(x: union_annot1, y: union_annot2) -> union_annot3:
            if isinstance(y, str):  # value = time
                return x.value  # returns <float>
            else:
                return x / y  # returns Quantity[m / s]

        annots = get_type_hints(union_func, include_extras=True)
        assert annots["x"] == union_annot1
        assert annots["y"] == union_annot2
        assert annots["return"] == union_annot3

    def test_quantity_subclass_typing(self):
        """Test type hint creation from a Quantity subclasses."""

        class Length(u.SpecificTypeQuantity):
            _equivalent_unit = u.m

        annot = Length[u.km]

        assert get_origin(annot) is Annotated
        assert get_args(annot) == (Length, u.km)
