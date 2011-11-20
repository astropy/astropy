# Licensed under a 3-clause BSD style license - see LICENSE.rst  

def test_config_file():
    from ..configs import get_config
    
    
    cfg1 = get_config('astropy')
    assert cfg1.filename.endswith('astropy.cfg')
    
    
    cfg2 = get_config('astropy.config')
    assert cg2.depth==1
    assert cfg2.name=='config'
    
def test_pkg_finder():
    from ..configs import _find_current_module
    
    assert _find_current_module(0).__name__ == 'astropy.config.configs'
    assert _find_current_module(1).__name__ == 'astropy.config.tests.test_configs'