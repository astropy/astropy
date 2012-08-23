from ..ui import read
from ..ipac import Ipac

DATA = '''
|   a  |   b   |
| char | char  |
ABBBBBBABBBBBBBA
'''

def test_ipac_default():
    # default should be right
    table = read(DATA, Reader=Ipac)
    assert table['a'][0] == 'ABBBBBB'
    assert table['b'][0] == 'ABBBBBBB'


def test_ipac_strict():
    # default should be strict
    table = read(DATA, Reader=Ipac, definition='strict')
    assert table['a'][0] == 'BBBBBB'
    assert table['b'][0] == 'BBBBBBB'
    
    
def test_ipac_left():
    # default should be strict
    table = read(DATA, Reader=Ipac, definition='left')
    assert table['a'][0] == 'BBBBBBA'
    assert table['b'][0] == 'BBBBBBBA'
    
    
def test_ipac_right():
    # default should be strict
    table = read(DATA, Reader=Ipac, definition='right')
    assert table['a'][0] == 'ABBBBBB'
    assert table['b'][0] == 'ABBBBBBB'