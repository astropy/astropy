# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ....tests.helper import raises

# LOCAL
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


def test_parse_vowarning(recwarn):
    config = {'pedantic': True,
              'filename': 'foo.xml'}
    pos = (42, 64)
    field = tree.Field(
        None, name='c', datatype='char',
        config=config, pos=pos)
    c = converters.get_converter(field, config=config, pos=pos)
    w = recwarn.pop(exceptions.W47)

    parts = exceptions.parse_vowarning(str(w.message))

    match = {
        'number'       : 47,
        'is_exception' : False,
        'nchar'        : 64,
        'warning'      : 'W47',
        'is_something' : True,
        'message'      : 'Missing arraysize indicates length 1',
        'doc_url'      : 'vo/api_exceptions.html#w47',
        'nline'        : 42,
        'is_warning'   : True
        }
    assert parts == match
