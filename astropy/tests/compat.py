import numpy as np

# assert_allclose doesn't exist in Numpy 1.4, so recreate it here
if tuple([int(x) for x in np.__version__.split(".", 2)[:2]]) < (1, 5):
    def assert_allclose(actual, desired, rtol=1e-7, atol=0,
                        err_msg='', verbose=True):
        from numpy.testing.utils import assert_array_compare

        def compare(x, y):
            return np.allclose(x, y, rtol=rtol, atol=atol)
        actual, desired = np.asanyarray(actual), np.asanyarray(desired)
        header = 'Not equal to tolerance rtol=%g, atol=%g' % (rtol, atol)
        assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
                             verbose=verbose, header=header)
else:
    from numpy.testing.utils import assert_allclose
