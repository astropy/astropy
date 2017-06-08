# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from . import Header, Card

from ...table import Table, QTable, Column
from ...time import Time


# Global time reference frame keywords
TIME_KEYWORDS = ['TIMESYS', 'MJDREF', 'JDREF', 'DATEREF', 'TREFPOS', 'TREFDIR', 'TIMEUNIT', 'TIMEOFFS']

# Column-specific time keywords
COLUMN_TIME_KEYWORDS = ['TCTYP[0-9]+',
                        'TRPOS[0-9]+',
                        'TUNIT[0-9]+']


def is_time_keyword(keyword):
    for c in TIME_KEYWORDS:
        if re.match(c, keyword) is not None:
            return True
    return False

def replace_time_table(table):
    '''
    Replace Time column in a Table with a normal column containing 
    each element as a vector of two doubles jd1 and jd2.
    jd = jd1 + jd2 represents time in the Julian Date format with high-precision.
    '''
    hdr = Header([Card(keyword = 'TIMESYS', value = 'UTC', comment = 'Default time scale'), Card(keyword = 'JDREF', value = 0.0, comment = 'Time columns are jd = jd1 + jd2')])
    
    time_cols = table.columns.isinstance(Time)
    for col in time_cols:
            table[col.info.name] = Column(np.empty(col.shape + (2,)))
            table[col.info.name][...,0] = col.jd1
            table[col.info.name][...,1] = col.jd2
            
            n = table.get_pos(col.info.name)
            hdr.append(Card(keyword = 'TCTYP%d' %n, value = col.scale))
            hdr.append(Card(keyword = 'APyFM%d' %n, value = col.format))                

    return table