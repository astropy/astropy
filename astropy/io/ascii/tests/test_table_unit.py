import pytest
from astropy.table import Table

def test_table_creation():
    """Check if a simple Table is created correctly."""
    t = Table({'a': [1, 2], 'b': [3, 4]})
    
    # Assertions to verify the data
    assert len(t) == 2
    assert list(t.columns.keys()) == ['a', 'b']
    assert t['a'][0] == 1
    assert t['b'][1] == 4
