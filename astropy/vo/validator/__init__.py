# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from ... import config as _config
from ...utils.data import get_pkg_data_contents


class _Conf(_config.ConfigNamespace):
    conesearch_master_list = _config.ConfigItem(
        'http://vao.stsci.edu/directory/NVORegInt.asmx/VOTCapabilityPredOpt?'
        'predicate=1%3D1&capability=conesearch&VOTStyleOption=2',
        'Cone Search services master list for validation.',
        aliases=['astropy.vo.validator.validate.cs_mstr_list']
    )

    conesearch_urls = _config.ConfigItem(
        get_pkg_data_contents(
            os.path.join('data', 'conesearch_urls.txt')).split(),
        'Only check these Cone Search URLs.',
        'list',
        aliases=['astropy.vo.validator.validate.cs_urls'])

    noncritical_warnings = _config.ConfigItem(
        ['W03', 'W06', 'W07', 'W09', 'W10', 'W15', 'W17', 'W20', 'W21', 'W22',
         'W27', 'W28', 'W29', 'W41', 'W42', 'W48', 'W50'],
        'VO Table warning codes that are considered non-critical',
        'list',
        aliases=['astropy.vo.validator.validate.noncrit_warnings'])
conf = _Conf()
