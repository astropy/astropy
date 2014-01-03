# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Defines constants used in `astropy.vo.samp`."""

__all__ = ['SAMP_STATUS_OK', 'SAMP_STATUS_WARNING', 'SAMP_STATUS_ERROR',
           'SAMP_HUB_SINGLE_INSTANCE', 'SAMP_HUB_MULTIPLE_INSTANCE',
           'SAMP_RESTRICT_GROUP', 'SAMP_RESTRICT_OWNER',
           'SAFE_MTYPES', 'SAMP_ICON']

__profile_version__ = "1.3"

#: General constant for samp.ok status string
SAMP_STATUS_OK = "samp.ok"
#: General constant for samp.warning status string
SAMP_STATUS_WARNING = "samp.warning"
#: General constant for samp.error status string
SAMP_STATUS_ERROR = "samp.error"

#: General constant to specify single instance Hub running mode
SAMP_HUB_SINGLE_INSTANCE = "single"
#: General constant to specify multiple instance Hub running mode
SAMP_HUB_MULTIPLE_INSTANCE = "multiple"

#: General constant to specify the access restriction (through Basic Authentication) to the GROUP
SAMP_RESTRICT_GROUP = "GROUP"
#: General constant to specify the access restriction (through Basic Authentication) to the OWNER
SAMP_RESTRICT_OWNER = "OWNER"

SAFE_MTYPES = ["samp.app.*", "samp.msg.progress", "table.*", "image.*",
               "coord.*", "spectrum.*", "bibcode.*", "voresource.*"]

SAMP_ICON = b"""
iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAMAAAAoLQ9TAAAABGdBTUEAALGPC/xhBQAAACBjSFJN
AAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAC2VBMVEX/gQD/hQD/hAD/kwD/
egD+ggD9eQD9eAD/jwD/mgD9aQD/iAH+gQL/jQD+dQD9bwD9LgD+MwD7AAD9LgT8JQD9KQn+PQD8
MRD8GwX+LgD+LQD9JAD5AAD9EwD9IgD9IQD9JQD/gwD/iQD/mAD/lwD/lwD/kAD/mgj/pBf+fgD+
hQL+gAD+jgD+jwD+kQH+kQD+qzr/hwD+jQD+kAD+hwD+fQH/hQD/iQP/hgD/hQD/jAr/qk7+zpf+
slv/hgD+hQL+ewH+fQH+gAL+dgD+dgD+ggn+gAb+dAH+bwD////+vH/+gAX+cwD+dQL+olL+/vv+
2rv+9Ov9fhb+dgX+bQL+u4r+x6L+9/D+sHn+bwf+bRL9WAD+jUv+y67+17/+ZAX+XAX+28j+ya7+
uJj+4M/+WgP9TAD++fX+iVj+0L/+ybD+Uwj9QgD++vj+f1L+VSD9OAD+VBT+PgL+4dj+q5X+UiD+
0sX+OwD+QA/+kHf++vf9Mwv+2tD+0MX9LQD+PQf9JwD+4t3+08f9IQf9hWz+tKj+qpz9aUv9IQD9
JAH9LwP9LQL8Mwz8GAD+5OD++PX9XEf8BQD8EAH8DQD8IAL8JgX9LAP9IgT9MhD8DAD9GQr9RTT9
HAb9GwX9GwP9JwH9JwT8AAD8AAD8CQj+7+79bWT8AAD9JBf9HAD9IwD+kQD+jwv/iQH/iwj+1Kz+
dAD9bQD+6dX+/fn+2LT+fxr+exn9cQD9ZAD////+jkD+lk7+5dT+9e3+dRv9XwD+kE79VwD+mlf+
j0f9XAD+49T+xqL9SQD9VgD9TwD9QwD+0rv+z7b+wqH+4ND9QgD9TgD+UgT+Xxj++Pb+y7T+/Pv+
u5z9NQD9QQD+Sgb+VBT+YTD+eEn+zrz+uJz9Qgr+glf9NAD9JAD+2s/+18n9NAL9IgD+6uH9NgD+
RCD+6+X/+/r+zsH+OxL9Kwb9Iwf9Nh/8EgD+mpH8DgCOXzztAAAAqnRSTlMAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAACARZVgIhsKgEBAxiU7v7++qoaAQNI7f79+/v8/eEfA1D9/vr9
/v7+/vi2JfH8/v7+/UOo+v7+/Igt9v7+/Kh0/P7+/KCg/P7+/HSo/P7+9i2I/P7++qhD/f7+/vzx
Jbb4/v7+/v36/v1QAx/h/fz7+/3+7UgDARqq+/7ulRgDAQEqbZODVRYBAgcNDmUAAAABYktHREmH
BeR8AAABG0lEQVQY0wEQAe/+AAABAiEiAyMkJSYnKAQpKgUABgcrCCwtLqovMKsxMjMJNAAKNQs2
N6w4OTo7PD2tPj8BAEAMQUJDREVGR0iuSUqvSw0ADkxNsE6xT7KztLVQUVK2UwAPVLdVuFa5uru8
vbdXWL5ZAFpbv0lcwMHCw8TFxl1ex18AYMhhYsnKuLjLzM3OY2TPZQBm0Gdo0dLTuNS41dZpatdr
AGzYbW7ZuNrb3N3eb0nfcHEAcuBzdOHi4+Tl5XXmdud3EAB44Xl6e+jp6ut8uH3sfn8RABKA7YGC
7oOEhYaHiImKE4sAFIyN746PkJGSk/CUlRWWFgCXF5iZmvFJm/KcnZ4YnxkaABugoRyio6Slpqcd
qKkeHyCF2XqFLvdwRQAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAxMy0xMi0xNlQxMToxNzo0NSswMTow
MIvg1+0AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMTMtMTItMTFUMjA6MzY6MzUrMDE6MDAb2LCcAAAA
AElFTkSuQmCC
"""

# TODO: document this global variable. Is this the right place for it?
_THREAD_STARTED_COUNT = 0
