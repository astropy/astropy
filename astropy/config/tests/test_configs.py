# Licensed under a 3-clause BSD style license - see LICENSE.rst  
    
def test_paths():
    from ..paths import get_config_dir, get_cache_dir
    
    assert 'astropy' in get_config_dir()
    assert 'astropy' in get_cache_dir()
    
def test_config_file():
    from ..configs import get_config,reload_config,save_config
    
    apycfg = get_config('astropy')
    assert apycfg.filename.endswith('astropy.cfg')
    
    cfgsec = get_config('astropy.config')
    assert cfgsec.depth==1
    assert cfgsec.name=='config'
    assert cfgsec.parent.filename.endswith('astropy.cfg')
    
    reload_config('astropy')
    #saving shouldn't change the file, because reload should have made sure it
    #is based on the current file
    save_config('astropy')
    
def test_pkg_finder():
    from ..configs import _find_current_module
    
    assert _find_current_module(0).__name__ == 'astropy.config.configs'
    assert _find_current_module(1).__name__ == 'astropy.config.tests.test_configs'