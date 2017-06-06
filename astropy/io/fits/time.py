# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

from . import Header, Card

from ...table import Table, QTable, Column
from ...time import Time


# Global time reference frame keywords
TIME_KEYWORDS = ['TIMESYS', 'MJDREF', 'JDREF', 'DATEREF', 'TREFPOS', 'TIMEUNIT', 'TIMEOFFS']

# Column-specific time keywords
COLUMN_TIME_KEYWORDS = ['TCTYP[0-9]+']

reference_date = dict('cxcsec' = ('1998-01-01 00:00:00','tt'), 'unix' = ('1970-01-01 00:00:00', 'utc'), 'gps' = ('1980-01-06 00:00:19', 'tai'))


def is_time_keyword(keyword):
    for c in TIME_KEYWORDS:
        if re.match(c, keyword) is not None:
            return True
    return False

def replace_time_table(table):
    
    hdr = Header()
    
    time_cols = table.columns.isinstance(Time)
    for col in time_cols:
        if col.format is in ('cxcsec', 'unix', 'gps'):
            if ref is not None:
                curr_scale = col.scale
                ref = Time(reference_date[col.format][0], scale = reference_date[col.format][1]).curr_scale
                c = Card(keyword = 'MJDREF', value = ref.mjd)
                hdr.append(c)
                continue
            
        else:
            table[col.info.name] = Column(np.empty(col.shape + (2,)))
            table[col.info.name][...,0] = col.jd1
            table[col.info.name][...,1] = col.jd2
            n = col.get_pos()
            c = Card(keyword = 'TCTYP%d' %n, value = col.scale)
            if col.format == 'jd':
                c = Card(keyword = 'JDREF', value = 0.0)
            elif col.format == 'mjd':
                c = Card(keyword = 'MJDREF', vaue = 0.0)
                

    return table