Configuration system Documentation
==================================

The astropy configuration system is designed to give users control of various
parameters used in astropy or affiliated packages without delving into the 
source code to make those changes.

.. todo::
    Move this portion of the documentation to a more appropriate place when a
    doc reorg happens

For Users
---------

To see the configuration options, look for your astropy configuration file.  
You can find it by doing::

    from astropy.config import get_config_dir
    
    print get_config_dir()

And you should see the location of your configuration directory. The standard 
scheme generally puts your configuration directory in ``$HOME/.astropy/config``,
but if you've set the environment variable `XDG_CONFIG_HOME` and the 
``$XDG_CONFIG_HOME/astropy`` directory exists, it will instead be there.

Once you've found the configuration file, open it with your favorite editor.
It should have all of the sections you might want, with descriptions and the 
type of the value that is accepted.  Feel free to edit this as you wish, and
any of these changes will be reflected when you next start Astropy.  Or, if you
want to see your changes immediately in your current Astropy session, just do::

    from astropy.config import reload_config()
    
    reload_config()

.. warning::
    The above is not true yet, because the setup doesn't automatically 
    populate the configuration files.  Hopefully it will be true soon, though.
    The :func:`astropy.config.configs._generate_all_config_items` function will
    already do this, basically, but there has to be some thought about how to
    make driver scripts that actually do this for each user, and coordinate
    when they get run so that everything is already built.

.. note::
    If for whatever reason your ``$HOME/.astropy`` directory is not accessible
    (i.e., you have astropy running somehow as root but you are not the root
    user), the best solution is to set the `XDG_CONFIG_HOME` and 
    `XDG_CACHE_HOME` environment variables pointing to directories, and create
    an ``astropy`` directory inside each of those.  Both the configuration and
    data download systems will then use those directories and never try to
    access the ``$HOME/.astropy`` directory.



For Developers
--------------
The Reference guide below describes the full interface - this is a summary of
typical practices.  The most common way to use the configuration system is as
follows::

    """ This is the docstring at the beginning of a module
    """
    from astropy.config import ConfigurationItem
    
    SOME_OPTION = ConfigurationItem('some_opt',1,'A description.')
    ANOTHER_OPTION = ConfigurationItem('anno_opt','a string val',
                            'A longer description of what this does.')  
                
    ... implementation ...
    def some_func():
        #to get the value of these options, I might do:
        something = SOME_OPTION()+2
        return ANOTHER_OPTION()+' Also, I added text.'
    
It is highly recommended that any configuration items be placed at the top of a
module like this, as they can then be easily found when viewing the source code
and the automated tools to generate the default configuration files can also
locate these items.



Reference/API
-------------
Below are the reference documentation for the three main config sub-packages:
`astropy.config.paths`, `astropy.config.configs`, and `astropy.config.data`. 
Note that all the public classes and functions in these sub-packages are 
imported into `astropy.config`, so the recommended usage is simply
``from astropy.config import foo``.

`astropy.config.paths`
^^^^^^^^^^^^^^^^^^^^^^

.. automodule::  astropy.config.paths
    :members:

`astropy.config.configs`
^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule::  astropy.config.configs
    :members:


`astropy.config.data`
^^^^^^^^^^^^^^^^^^^^^

.. automodule::  astropy.config.data
    :members: