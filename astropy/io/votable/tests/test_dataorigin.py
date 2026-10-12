# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

import astropy.io.votable
import astropy.io.votable.dataorigin as dataorigin
from astropy.table import Column, Table
from astropy.utils.exceptions import AstropyDeprecationWarning


def __generate_votable_test():
    table = Table(
        [
            Column(name="id", data=[1, 2, 3, 4]),
            Column(name="bmag", unit="mag", data=[5.6, 7.9, 12.4, 11.3]),
        ]
    )
    return astropy.io.votable.tree.VOTableFile().from_table(table)


__TEST_PUBLISHER_NAME = "astropy"
__TEST_CREATOR1_NAME = "her"
__TEST_CREATOR2_NAME = "him"
__TEST_CONTACT = "me@mail.fr"
__TEST_SERVICE_PROTOCOL = "ivo://std"
__TEST_DATA_IVOID = "ivo://id"
__TEST_JOURNAL = "AJ"


def __add_origin():
    vot = __generate_votable_test()
    dataorigin.add_data_origin_info(vot, "publisher", __TEST_PUBLISHER_NAME)
    dataorigin.add_data_origin_info(
        vot, "contact", __TEST_CONTACT, "Contact email address"
    )
    dataorigin.add_data_origin_info(
        vot.resources[0], "creator", __TEST_CREATOR1_NAME, "Author name"
    )
    dataorigin.add_data_origin_info(
        vot.resources[0], "creator", __TEST_CREATOR2_NAME, "Author name"
    )
    dataorigin.add_data_origin_info(
        vot.resources[0], "data_ivoid", __TEST_DATA_IVOID, "ivoid"
    )
    dataorigin.add_data_origin_info(
        vot.resources[0], "journal", __TEST_JOURNAL, "journal abbrev"
    )

    return vot


def test_data_origin_update():
    vot = __add_origin()
    do = dataorigin.extract_data_origin(vot)
    do.update_votable()

    vot = __generate_votable_test()
    dataorigin.add_data_origin_info(vot, "publisher", __TEST_PUBLISHER_NAME, "Pub.")
    do = dataorigin.extract_data_origin(vot)
    # do.query.publisher = __TEST_PUBLISHER_NAME
    do.query.standardID = __TEST_SERVICE_PROTOCOL
    do.origin = [dataorigin.DatasetOrigin(vot.resources[0])]
    do.origin[0].creator = __TEST_CREATOR1_NAME

    do.update_votable()
    assert len(vot.infos) == 2
    assert len(vot.resources[0].infos) == 1
    assert do.query.service_protocol == __TEST_SERVICE_PROTOCOL


def test_dataorigin():
    vot = __add_origin()
    do = dataorigin.extract_data_origin(vot)

    assert do.query.publisher == __TEST_PUBLISHER_NAME
    assert len(do.origin) == 1 and len(do.origin[0].creator) == 2
    assert do.origin[0].creator[0] == __TEST_CREATOR1_NAME
    assert len(str(do)) > 1
    assert (do.origin[0].get_votable_element()) is not None
    assert do.origin[0].journal[0] == __TEST_JOURNAL
    assert do.origin[0].data_ivoid[0] == __TEST_DATA_IVOID
    # the deprecated names stay backward compatible and redirect to the new ones
    with pytest.warns(AstropyDeprecationWarning, match="data_ivoid"):
        assert do.origin[0].ivoid[0] == __TEST_DATA_IVOID
    with pytest.warns(AstropyDeprecationWarning, match="journal"):
        assert do.origin[0].editor[0] == __TEST_JOURNAL

    dores = dataorigin.extract_data_origin(vot.resources[0])
    dot = dataorigin.extract_data_origin(vot.resources[0].tables[0])

    assert dores.origin[0].creator == do.origin[0].creator
    assert len(dot.origin) == 0

    with pytest.raises(ValueError, match="QueryOrigin contact already exists"):
        # QueryOrigin attributes are unique
        dataorigin.add_data_origin_info(vot, "contact", "othercontact")

    with pytest.raises(TypeError, match="Bad type of vot_element"):
        # contact should be global in VOFile node
        dataorigin.add_data_origin_info(vot.resources[0], "contact", "othercontact")

    with pytest.raises(ValueError, match="Unknown DataOrigin info name."):
        dataorigin.add_data_origin_info(vot.resources[0], "unexistingkey", "value")

    for resource in do:
        assert len(resource.creator) > 0


def test_dataorigin_deprecated_attributes():
    # The 1.2 names ``ivoid`` and ``editor`` were renamed to ``data_ivoid`` and
    # ``journal``. They must still work for backward compatibility while emitting
    # a deprecation warning that points at the new name.
    origin = dataorigin.DatasetOrigin()
    origin.data_ivoid = ["ivo://id"]
    origin.journal = ["AJ"]

    with pytest.warns(AstropyDeprecationWarning, match="data_ivoid"):
        assert origin.ivoid == ["ivo://id"]
    with pytest.warns(AstropyDeprecationWarning, match="journal"):
        assert origin.editor == ["AJ"]

    # the deprecated setters must redirect to the new attributes
    with pytest.warns(AstropyDeprecationWarning, match="data_ivoid"):
        origin.ivoid = ["ivo://other"]
    assert origin.data_ivoid == ["ivo://other"]

    with pytest.warns(AstropyDeprecationWarning, match="journal"):
        origin.editor = ["ApJ"]
    assert origin.journal == ["ApJ"]


def test_add_data_origin_info_deprecated_names():
    # Passing an obsolete INFO name to add_data_origin_info must emit a
    # deprecation warning and store the value under the current name.
    vot = __generate_votable_test()
    with pytest.warns(AstropyDeprecationWarning, match="data_ivoid"):
        dataorigin.add_data_origin_info(vot.resources[0], "ivoid", __TEST_DATA_IVOID)
    with pytest.warns(AstropyDeprecationWarning, match="journal"):
        dataorigin.add_data_origin_info(vot.resources[0], "editor", __TEST_JOURNAL)

    # the obsolete names are translated to the current vocabulary
    names = [info.name for info in vot.resources[0].infos]
    assert "data_ivoid" in names and "journal" in names
    assert "ivoid" not in names and "editor" not in names

    do = dataorigin.extract_data_origin(vot)
    assert do.origin[0].data_ivoid[0] == __TEST_DATA_IVOID
    assert do.origin[0].journal[0] == __TEST_JOURNAL


def test_dataorigin_unsupported_input_error():
    table = Table.read(__generate_votable_test())

    with pytest.raises(TypeError, match="input vot_element type is not supported."):
        dataorigin.extract_data_origin(table)
