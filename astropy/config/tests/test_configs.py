def test_config_file():
    from ..configs import get_config
    
    
    cfg = get_config('config')
    assert cfg.filename.endswith('config.cfg')