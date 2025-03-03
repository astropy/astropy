# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.io.votable
import astropy.io.votable.dataorigin as dataorigin
from astropy.table import Table, Column


def __generate_votable_test():
    table = Table([Column(name='id', data=[1, 2, 3, 4]),
                   Column(name='bmag', unit='mag', data=[5.6, 7.9, 12.4, 11.3])])
    return astropy.io.votable.tree.VOTableFile().from_table(table)


__TEST_PUBLISHER_NAME = "astropy"
__TEST_CREATOR_NAME = "me"
__TEST_CONTACT = "me@mail.fr"


def __add_origin():
    vot = __generate_votable_test()
    dataorigin.add_data_origin_info(vot, "publisher", __TEST_PUBLISHER_NAME)
    dataorigin.add_data_origin_info(vot, "contact", __TEST_CONTACT, "Contact email address")
    dataorigin.add_data_origin_info(vot.resources[0], "creator", __TEST_CREATOR_NAME, "Author name")
    return vot


def test_dataorigin():
    vot = __add_origin()
    do = dataorigin.extract_data_origin(vot)
    assert do.query.publisher == __TEST_PUBLISHER_NAME
    assert len(do.origin) == 1 and len(do.origin[0].creator) == 1
    assert do.origin[0].creator[0] == __TEST_CREATOR_NAME

