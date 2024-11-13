# Licensed under a 3-clause BSD style license - see LICENSE.rst
# LOCAL
import io

from astropy.io.votable import parse
from astropy.utils.data import get_pkg_data_filename


def test_resource_groups():
    # Read the VOTABLE
    votable = parse(get_pkg_data_filename("data/resource_groups.xml"))

    resource = votable.resources[0]
    groups = resource.groups
    params = resource.params

    # Test that params inside groups are not outside

    assert len(groups[0].entries) == 1
    assert groups[0].entries[0].name == "ID"

    assert len(params) == 2
    assert params[0].name == "standardID"
    assert params[1].name == "accessURL"


def test_roundtrip():
    # Issue #16511 VOTable writer does not write out GROUPs within RESOURCEs

    # Read the VOTABLE
    votable = parse(get_pkg_data_filename("data/resource_groups.xml"))

    bio = io.BytesIO()
    votable.to_xml(bio)
    bio.seek(0)
    votable = parse(bio)

    resource = votable.resources[0]
    groups = resource.groups
    params = resource.params

    # Test that params inside groups are not outside

    assert len(groups[0].entries) == 1
    assert groups[0].entries[0].name == "ID"

    assert len(params) == 2
    assert params[0].name == "standardID"
    assert params[1].name == "accessURL"
