# Licensed under a 3-clause BSD style license - see LICENSE.rst  

def test_config_file():
    from ..configs import get_config
    
    
    cfg = get_config('config')
    assert cfg.filename.endswith('config.cfg')