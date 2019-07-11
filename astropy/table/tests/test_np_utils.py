import numpy as np

from astropy.table import np_utils


def test_common_dtype():
    """
    Test that allowed combinations are those expected.
    """
    dtype = [('int', int),
             ('uint8', np.uint8),
             ('float32', np.float32),
             ('float64', np.float64),
             ('str', 'S2'),
             ('uni', 'U2'),
             ('bool', bool),
             ('object', np.object_)]
    arr = np.empty(1, dtype=dtype)
    fail = set()
    succeed = set()
    for name1, type1 in dtype:
        for name2, type2 in dtype:
            try:
                np_utils.common_dtype([arr[name1], arr[name2]])
                succeed.add(f'{name1} {name2}')
            except np_utils.TableMergeError:
                fail.add(f'{name1} {name2}')

    # known bad combinations
    bad = set(['str int', 'str bool', 'uint8 bool', 'uint8 str', 'object float32',
               'bool object', 'uni uint8', 'int str', 'bool str', 'bool float64',
               'bool uni', 'str float32', 'uni float64', 'uni object', 'bool uint8',
               'object float64', 'float32 bool', 'str uint8', 'uni bool', 'float64 bool',
               'float64 object', 'int bool', 'uni int', 'uint8 object', 'int uni', 'uint8 uni',
               'float32 uni', 'object uni', 'bool float32', 'uni float32', 'object str',
               'int object', 'str float64', 'object int', 'float64 uni', 'bool int',
               'object bool', 'object uint8', 'float32 object', 'str object', 'float64 str',
               'float32 str'])
    assert fail == bad

    good = set(['float64 int', 'int int', 'uint8 float64', 'uint8 int', 'str uni',
                'float32 float32', 'float64 float64', 'float64 uint8', 'float64 float32',
                'int uint8', 'int float32', 'uni str', 'int float64', 'uint8 float32',
                'float32 int', 'float32 uint8', 'bool bool', 'uint8 uint8', 'str str',
                'float32 float64', 'object object', 'uni uni'])
    assert succeed == good
