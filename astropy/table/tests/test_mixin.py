from ..table import QTable
from ... import units as u

# ISSUES
# - Overlap with SkyCoord.name => SkyCoord.frame.name => 'icrs' by default.
#   This gets clobbered by directly setting `frame` attribute.

def test_attributes(mixin_cols):
    """
    Required attributes for a column can be set.
    """
    m = mixin_cols['m']
    m.name = 'a'
    assert m.name == 'a'

    m.description = 'a'
    assert m.description == 'a'

    if not isinstance(m, u.Quantity):
        m.unit = u.m
    assert m.unit is u.m

    if hasattr(m, '_table_format'):
        m._table_format = 'a'
        assert m._table_format == 'a'
    else:
        m.format = 'a'
        assert m.format == 'a'

    m.meta = {'a': 1}
    assert m.meta == {'a': 1}


def test_make_table(table_types, mixin_cols):
    """
    Make a table with the columns in mixin_cols, which is an ordered dict of
    three cols: 'a' and 'b' are table_types.Column type, and 'm' is a mixin.
    """
    t = table_types.Table(mixin_cols)
    m = mixin_cols['m']
    tm = t['m']
    if isinstance(m, u.Quantity) and type(t) is not QTable:
        assert type(tm) is table_types.Column
    else:
        assert type(tm) is type(m)
