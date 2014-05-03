from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from ...tests.helper import catch_warnings
from ..data import get_pkg_data_filename
from ...config import configuration


def test_alias():
    from ...cosmology import core, WMAP9, WMAP7

    # REMOVE in astropy 0.5

    with catch_warnings() as w:
        x = core.DEFAULT_COSMOLOGY()
    assert x == WMAP9
    assert len(w) == 1
    assert str(w[0].message) == (
        "'astropy.cosmology.core.DEFAULT_COSMOLOGY' is deprecated, "
        "and is no longer defined as a configuration item. Use "
        "'astropy.cosmology.core.default_cosmology.get()' instead.")

    with catch_warnings() as w:
        core.DEFAULT_COSMOLOGY.set('WMAP7')
    assert core.default_cosmology.get() == WMAP7
    assert len(w) == 1
    assert str(w[0].message) == (
        "'astropy.cosmology.core.DEFAULT_COSMOLOGY' is deprecated, "
        "and is no longer defined as a configuration item. Use "
        "'astropy.cosmology.core.default_cosmology.set()' instead.")

    with catch_warnings() as w:
        with core.DEFAULT_COSMOLOGY.set_temp('WMAP9'):
            assert core.default_cosmology.get() == WMAP9
        assert core.default_cosmology.get() == WMAP7
    assert len(w) == 1
    assert str(w[0].message) == (
        "'astropy.cosmology.core.DEFAULT_COSMOLOGY' is deprecated, "
        "and is no longer defined as a configuration item. Use "
        "'astropy.cosmology.core.default_cosmology.set_temp()' instead.")

    with catch_warnings() as w:
        core.DEFAULT_COSMOLOGY.reload()
    assert core.default_cosmology.get() == WMAP9
    assert len(w) == 1
    assert str(w[0].message) == (
        "'astropy.cosmology.core.DEFAULT_COSMOLOGY' is deprecated, "
        "and is no longer defined as a configuration item.")


class TestAliasRead(object):
    def setup_class(self):
        configuration._override_config_file = get_pkg_data_filename('data/alias.cfg')

    def test_alias_read(self):
        from ...cosmology import core, WMAP9, WMAP7

        with catch_warnings() as w:
            core.DEFAULT_COSMOLOGY.reload()
            assert core.default_cosmology.get() == WMAP7

        assert len(w) == 1
        assert str(w[0].message) == (
            "'astropy.cosmology.core.DEFAULT_COSMOLOGY' is deprecated, and is "
            "no longer defined as a configuration item.")

    def teardown_class(self):
        from astropy.utils.data import conf

        configuration._override_config_file = None
        conf.reload()
