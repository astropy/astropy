# Licensed under a 3-clause BSD style license - see LICENSE.rst

# LOCAL
from ....tests.helper import catch_warnings

from .. import converters
from .. import exceptions
from .. import tree


def test_reraise():
    def fail():
        raise RuntimeError("This failed")

    try:
        try:
            fail()
        except RuntimeError as e:
            exceptions.vo_reraise(e, additional="From here")
    except RuntimeError as e:
        assert "From here" in str(e)
    else:
        assert False


def test_parse_vowarning():
    config = {'pedantic': True,
              'filename': 'foo.xml'}
    pos = (42, 64)
    with catch_warnings(exceptions.W47) as w:
        field = tree.Field(
            None, name='c', datatype='char',
            config=config, pos=pos)
        c = converters.get_converter(field, config=config, pos=pos)

    parts = exceptions.parse_vowarning(str(w[0].message))

    match = {
        'number': 47,
        'is_exception': False,
        'nchar': 64,
        'warning': 'W47',
        'is_something': True,
        'message': 'Missing arraysize indicates length 1',
        'doc_url': 'io/votable/api_exceptions.html#w47',
        'nline': 42,
        'is_warning': True
        }
    assert parts == match
