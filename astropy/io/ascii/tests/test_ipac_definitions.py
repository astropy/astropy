# Licensed under a 3-clause BSD style license - see LICENSE.rst
from ..ui import read
from ..ipac import Ipac

DATA = '''
|   a  |   b   |
| char | char  |
ABBBBBBABBBBBBBA
'''

def test_ipac_default():
    # default should be ignore
    table = read(DATA, Reader=Ipac)
    assert table['a'][0] == 'BBBBBB'
    assert table['b'][0] == 'BBBBBBB'


def test_ipac_ignore():
    table = read(DATA, Reader=Ipac, definition='ignore')
    assert table['a'][0] == 'BBBBBB'
    assert table['b'][0] == 'BBBBBBB'
    
    
def test_ipac_left():
    table = read(DATA, Reader=Ipac, definition='left')
    assert table['a'][0] == 'BBBBBBA'
    assert table['b'][0] == 'BBBBBBBA'
    
    
def test_ipac_right():
    table = read(DATA, Reader=Ipac, definition='right')
    assert table['a'][0] == 'ABBBBBB'
    assert table['b'][0] == 'ABBBBBBB'