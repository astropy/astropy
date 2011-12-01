# Licensed under a 3-clause BSD style license - see LICENSE.rst  

def test_paths():
    from ..paths import get_config_dir, get_cache_dir
    
    assert 'astropy' in get_config_dir()
    assert 'astropy' in get_cache_dir()
    
def test_config_file():
    from ..configs import get_config,reload_config,save_config
    from os.path import exists
    
    apycfg = get_config('astropy')
    assert apycfg.filename.endswith('astropy.cfg')
    
    cfgsec = get_config('astropy.config')
    assert cfgsec.depth==1
    assert cfgsec.name=='config'
    assert cfgsec.parent.filename.endswith('astropy.cfg')
    
    reload_config('astropy')
    
    #saving shouldn't change the file, because reload should have made sure it
    #is based on the current file.  But don't do it if there's no file
    if exists(apycfg.filename):
        save_config('astropy')

def test_configitem(tmpdir):
    from ..configs import ConfigurationItem,get_config
    from shutil import copy
    
    ci = ConfigurationItem('tstnm',34,'this is a Description')

    assert ci.module=='astropy.config.tests.test_configs'
    assert ci()==34
    assert ci.description=='this is a Description'
    
    sec = get_config(ci.module)
    assert sec['tstnm'] == 34
    assert sec.comments['tstnm'][0] == 'this is a Description'
    assert sec.comments['tstnm'][1] == 'integer'
    
    ci.description = 'updated Descr'
    ci.set(32)
    assert ci()==32
    assert sec.comments['tstnm'][0] == 'updated Descr'
    
    #now try saving
    apycfg = sec
    while apycfg.parent is not apycfg:
        apycfg = apycfg.parent
    f = tmpdir.join('astropy.cfg')
    apycfg.write(f)
    lns = f.read().split('\n')
    assert 'tstnm = 32' in lns
    assert '# updated Descr' in lns
    
    oldfn = apycfg.filename
    try:
        tmpdir.join('astropy.cfg').copy(tmpdir.join('astropy.cfg2'))
        apycfg.filename = str(tmpdir.join('astropy.cfg2').realpath())
        
        ci.set(30)
        ci.save()
        
        with open(str(apycfg.filename)) as f:
            lns = f.readlines()
            assert '[config.tests.test_configs]\n' in lns
            assert 'tstnm = 30\n' in lns
        
        ci.save(31)
        
        with open(str(apycfg.filename)) as f:
            lns = f.readlines()
            assert '[config.tests.test_configs]\n' in lns
            assert 'tstnm = 31\n' in lns
            
        #also try to save one that doesn't yet exist
        apycfg.filename = str(tmpdir.join('astropy.cfg3').realpath())
        ci.save()
        
        with open(str(apycfg.filename)) as f:
            lns = f.readlines()
            assert '[config.tests.test_configs]\n' in lns
            assert 'tstnm = 30\n' in lns
            
    finally:
        apycfg.filename = oldfn
        
    #also try to save one that doesn't yet exist
    
def test_configitem_types():
    from ..configs import ConfigurationItem
    from pytest import raises
    
    ci1 = ConfigurationItem('tstnm1',34)
    assert isinstance(ci1(),int)
    
    ci2 = ConfigurationItem('tstnm2',34.3)
    assert isinstance(ci2(),float)
    
    ci3 = ConfigurationItem('tstnm3',True)
    assert isinstance(ci3(),bool)
    
    ci4 = ConfigurationItem('tstnm4','astring')
    assert isinstance(ci4(),str)
    
    with raises(TypeError):
        ci1.set(34.3)
    ci2.set(12) #this would should succeed as up-casting
    with raises(TypeError):
        ci3.set('fasd')
    with raises(TypeError):
        ci4.set(546.245)
    
    
def test_configitem_options(tmpdir):
    from ..configs import ConfigurationItem,get_config
    from pytest import raises
    
    cio = ConfigurationItem('tstnmo',['op1','op2','op3'])
    sec = get_config(cio.module)
    
    assert isinstance(cio(),str)
    assert cio() == 'op1'
    assert sec['tstnmo'] == 'op1'
    
    
    cio.set('op2')
    with raises(TypeError):
        cio.set('op5')
    assert sec['tstnmo'] == 'op2'
    
    #now try saving
    apycfg = sec
    while apycfg.parent is not apycfg:
        apycfg = apycfg.parent
    f = tmpdir.join('astropy.cfg')
    apycfg.write(f)
    lns = f.read().split('\n')
    
    assert '# option(op1, op2, op3)' in lns
    assert 'tstnmo = op2' in lns

    