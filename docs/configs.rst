Configuration system Documentation (`astropy.config`)
=====================================================

Introduction
------------

The astropy configuration system is designed to give users control of various
parameters used in astropy or affiliated packages without delving into the
source code to make those changes.


Getting Started
---------------

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
    populate the configuration files. Hopefully it will be true soon, though.
    The :func:`~astropy.config.configuration._generate_all_config_items`
    function will already do this, basically, but there has to be some
    thought about how to make driver scripts that actually do this for
    each user, and coordinate when they get run so that everything is
    already built.

.. note::
    If for whatever reason your ``$HOME/.astropy`` directory is not accessible
    (i.e., you have astropy running somehow as root but you are not the root
    user), the best solution is to set the `XDG_CONFIG_HOME` and
    `XDG_CACHE_HOME` environment variables pointing to directories, and create
    an ``astropy`` directory inside each of those.  Both the configuration and
    data download systems will then use those directories and never try to
    access the ``$HOME/.astropy`` directory.


Advanced usage of the configuration system
------------------------------------------

Changing Values at Run-time
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The configuration system is most conviniently used, by modifying
configuration files as described above. Values can also, however, be
modified in an active python session using the
:meth:`~astropy.config.configuration.ConfigurationItem.set` method. A run-time
`ConfigurationItem` object can be used to make these changes. These items
are found in the same module as the configuration section they are in,
and usually have the same name as in the configuration files, but in all
caps. Alternatively, they may be located with the
:func:`~astropy.config.configuration.get_config_items` function.

For example, if there is a part of your configuration file that looks like::

    [config.data]
    # URL for astropy remote data site.
    dataurl = http://data.astropy.org/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 3.0


You should be able to modify the values at run-time this way::

    from astropy.config.data import DATAURL, REMOTE_TIMEOUT

    DATAURL.set('http://astropydata.mywebsite.com')
    REMOTE_TIMEOUT.set(4.5)

Or alternatively::

    from astropy.config import get_config_items

    items = get_config_items('astropy.config.data')
    items['DATAURL'].set('http://astropydata.mywebsite.com')
    items['REMOTE_TIMEOUT'].set('4.5')

Note that this will *not* permanently change these values in the configuration
files - just for the current session.  To modify the values in the configuration
files, you can do::

    DATAURL.save()
    REMOTE_TIMEOUT.save()

Or to save all modifications to configuration items in `astropy.config.data`
(which includes the changes made above), do::

    from astropy.config import save_config

    save_config('astropy.config.data')



Developer Usage
^^^^^^^^^^^^^^^

The most common way to use the configuration system in a module is as follows::

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

It is highly recommended that any configuration items be placed at the
top of a module like this, as they can then be easily found when viewing
the source code. Additionally, the automated tools to generate the
default configuration files can also locate these items, and will use the
description here to produce descriptive comments in the configuration file.

It is important to realize that whil the typical usage of the configuration
system uses the astropy configuration files, the configuration items can also
be changed at runtime.  Thus, in general, you should access the value of an item
via the ``SOME_OPTION()`` call syntax above.  This always accesses an up-to-date
value for the item rather than assuming that its value remains fixed until the
module is next imported.  If you need to use an item as an initialization
parameter (such that changing it at runtime has no effect), be sure to note this
in an item's description.


See Also
--------

:doc:`logging` (overview of `astropy.config.logging`)


API/Reference
-------------

.. automodapi:: astropy.config
    :no-inheritance-diagram:






