# Licensed under a 3-clause BSD style license - see LICENSE.rst

from functools import wraps

import pytest

from astropy.utils.compat.optional_deps import HAS_PYTEST_MPL


def figure_test(*args, **kwargs):
    """
    A decorator that defines a figure test.

    This automatically decorates tests with mpl_image_compare with common
    options used by all figure tests in astropy, and also adds the decorator
    to allow remote data to be accessed.
    """
    # NOTE: the savefig_kwargs option below is to avoid using PNG files with
    # the matplotlib version embedded since this changes for every developer
    # version.

    tolerance = kwargs.pop("tolerance", 0)
    style = kwargs.pop("style", {})
    savefig_kwargs = kwargs.pop("savefig_kwargs", {})
    savefig_kwargs["metadata"] = {"Software": None}

    def decorator(test_function):
        @pytest.mark.remote_data
        @pytest.mark.mpl_image_compare(
            tolerance=tolerance, style=style, savefig_kwargs=savefig_kwargs, **kwargs
        )
        @pytest.mark.skipif(
            not HAS_PYTEST_MPL, reason="pytest-mpl is required for the figure tests"
        )
        @wraps(test_function)
        def test_wrapper(*args, **kwargs):
            return test_function(*args, **kwargs)

        return test_wrapper

    # If the decorator was used without any arguments, the only positional
    # argument will be the test to decorate so we do the following:
    if len(args) == 1:
        return decorator(*args)

    return decorator
