.. _astropy_config:

***************************************
Configuration system (`astropy.config`)
***************************************

Introduction
============

The astropy configuration system is designed to give users control of various
parameters used in astropy or affiliated packages without delving into the
source code to make those changes.


Getting Started
===============

The Astropy configuration options are most easily set by modifying the
configuration file.  It will be automatically generated with all the
default values the first time you import Astropy.  You can find the
exact location by doing::

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

    from astropy.config import reload_config

    reload_config()

.. note::
    If for whatever reason your ``$HOME/.astropy`` directory is not accessible
    (i.e., you have astropy running somehow as root but you are not the root
    user), the best solution is to set the `XDG_CONFIG_HOME` and
    `XDG_CACHE_HOME` environment variables pointing to directories, and create
    an ``astropy`` directory inside each of those.  Both the configuration and
    data download systems will then use those directories and never try to
    access the ``$HOME/.astropy`` directory.


Using `config`
==============

Accessing Values
----------------

The value of a configuration object can be accessed by calling the
object as a function.  For instance to get the default URL for
astropy remote data do::

    >>> from astropy.utils.data import DATAURL
    >>> DATAURL()
    'http://data.astropy.org/'

To interactively see more information about a configuration object there are
two options::

    >>> DATAURL
    <ConfigurationItem: name='dataurl' value='http://data.astropy.org/' at 0x1032c4290>
    >>> print DATAURL
    ConfigurationItem: dataurl
      cfgtype='string'
      defaultvalue='http://data.astropy.org/'
      description='URL for astropy remote data site.'
      module=astropy.utils.data
      value='http://data.astropy.org/'

Changing Values at Run-time
---------------------------

The configuration system is most conveniently used by modifying
configuration files as described above. Values can also, however, be
modified in an active python session using the
:meth:`~astropy.config.configuration.ConfigurationItem.set` method. A run-time
`ConfigurationItem` object can be used to make these changes. These items
are found in the same module as the configuration section they are in,
and usually have the same name as in the configuration files, but in all
caps. Alternatively, they may be located with the
:func:`~astropy.config.configuration.get_config_items` function.

For example, if there is a part of your configuration file that looks like::

    [utils.data]

    # URL for astropy remote data site.
    dataurl = http://data.astropy.org/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 3.0


You should be able to modify the values at run-time this way::

    >>> from astropy.utils.data import DATAURL, REMOTE_TIMEOUT
    >>> DATAURL()
    'http://data.astropy.org/'
    >>> DATAURL.set('http://astropydata.mywebsite.com')
    >>> DATAURL()
    'http://astropydata.mywebsite.com'
    >>> REMOTE_TIMEOUT()
    3.0
    >>> REMOTE_TIMEOUT.set(4.5)
    >>> REMOTE_TIMEOUT()
    4.5

Or alternatively::

    >>> from astropy.config import get_config

    >>> items = get_config('astropy.utils.data')
    >>> items['dataurl'].set('http://astropydata.mywebsite.com')
    >>> items['remote_timeout'].set('4.5')

Note that this will *not* permanently change these values in the configuration
files - just for the current session.  To change the configuration files,
after you've made your changes, you can do::

    >>> DATAURL.save()
    >>> REMOTE_TIMEOUT.save()

Or to save all modifications to configuration items in `astropy.utils.data`
(which includes the changes made above), do::

    >>> from astropy.config import save_config
    >>> save_config('astropy.utils.data')

Reloading Configuration
-----------------------

Instead of modifying the variables in python, you can also modify the
configuration files and then reload them.  For example, if you modify the
configuration file to say::

    [utils.data]

    # URL for astropy remote data site.
    dataurl = http://myotherdata.mywebsite.com/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 6.3

And then run the following commands::

    >>> DATAURL.reload()
    >>> REMOTE_TIMEOUT.reload()

This should update the variables with the values from the configuration file::

    >>> DATAURL()
    'http://myotherdata.mywebsite.com/'
    >>> REMOTE_TIMEOUT()
    6.3

Or if you want to reload all astropy configuration at once, use the
`~astropy.config.configuration.reload_config` function::

    >>> config.reload_config('astropy')


Developer Usage
---------------

Configuration items should be used wherever an option or setting is
needed that is either tied to a system configuration or should persist
across sessions of astropy or an affiliated package. Admittedly, this is
only a guideline, as the precise cases where a configuration item is
preferred over, say, a keyword option for a function is somewhat personal
preference. It is the preferred form of persistent configuration,
however, and astropy packages must all use it (and it is recommended for
affiliated packages).

The Reference guide below describes the full interface for a
`ConfigurationItem` - this is a guide for *typical* developer usage. In
almost all cases, a configuration item should be defined and used in the
following manner::

    """ This is the docstring at the beginning of a module
    """
    from astropy.config import ConfigurationItem

    SOME_SETTING = ConfigurationItem('some_setting', 1, 'A description.')
    ANOTHER_SETTING = ConfigurationItem('another_set', 'a string val',
                                       'A longer description of what this does.')

    ... implementation ...
    def some_func():
        #to get the value of these options, I might do:
        something = SOME_SETTING() + 2
        return ANOTHER_SETTING() + ' Also, I added text.'

It is highly recommended that any configuration items be placed at the
top of a module like this, as they can then be easily found when viewing
the source code and the automated tools to generate the default
configuration files can also locate these items.

Item Types and Validation
^^^^^^^^^^^^^^^^^^^^^^^^^

If not otherwise specified, a `ConfigurationItem` gets its type from the type of
the `defaultvalue` it is given when it is created.  The item can only be set
to be an object of this type.  Hence::

    SOME_SETTING = ConfigurationItem('some_setting', 1, 'A description.')
    SOME_SETTING.set(1.2)

will fail, because ``1.2`` is a float and ``1`` is an int.

Note that if you want the configuration item to be limited to a particular set
of specific options, you should pass in a list as the `defaultvalue` option.
The first entry in the list will be taken as the default, and the list as a whole
gives all the valid options.  For example::

    AN_OPTION = ConfigurationItem('an_option', ['a', 'b', 'c'], "This option can be 'a', 'b', or 'c'")
    AN_OPTION.set('b')  # succeeds
    AN_OPTION.set('c')  # succeeds
    AN_OPTION.set('d')  # fails!
    AN_OPTION.set(6)  # fails!

Finally, a ConfigurationItem can be explicitly give a type via the `cfgtype` option::

    AN_INT_SETTING = ConfigurationItem('an_int_setting', 1, 'A description.', cfgtype='integer')
    AN_INT_SETTING.set(3)  # works fine
    AN_INT_SETTING.set(4.2)  #fails!

If the default value's type doesn't match `cfgtype`, the `ConfigurationItem` cannot be created::

    >>> AN_INT_SETTING = ConfigurationItem('an_int_setting', 4.2, 'A description.', cfgtype='integer')
    VdtTypeError: the value "4.2" is of the wrong type.

hence the default behavior (of automatically determining `cfgtype`) is usually what
you want.  The main exception is when you want your configuration item to be a list.
The default behavior will treat that as a list of *options* unless you explicitly
tell it that the `ConfigurationItem` itself is supposed to be a list::

    >>> A_LIST_SETTING = ConfigurationItem(a_list_setting', [1, 2, 3], 'A description.')
    >>> A_LIST_SETTING()
    1
    >>> A_LIST_SETTING = ConfigurationItem(a_list_setting', [1, 2, 3], 'A description.', cfgtype='list')
    >>> A_LIST_SETTING()
    [1, 2, 3]

Details of all the valid `cfgtype` items can be found in the
`validation section of the configobj  manual <http://www.voidspace.org.uk/python/validate.html#the-standard-functions>`_.
We simply list the valid values here for quick reference:

* 'integer'
* 'float'
* 'boolean'
* 'string'
* 'ip_addr'
* 'list'
* 'tuple'
* 'int_list'
* 'float_list'
* 'bool_list'
* 'string_list'
* 'ip_addr_list'
* 'mixed_list'
* 'option'
* 'pass'




Usage Tips
^^^^^^^^^^

There are a couple important gotchas to remember about using configuration
items in your code. First, it is tempting to do something like::

    SOME_SETTING = ConfigurationItem('SOME_SETTING' ,1 ,'A description.')

    def some_func():
        return SOME_SETTING + 2  # WRONG, you wanted SOME_SETTING() + 2

but this is incorrect, because ``SOME_SETTING`` instead of
``SOME_SETTING()`` will yield a `ConfigurationItem` object, instead of the
*value* of that item (an integer, in this case).

The second point to keep in mind is that `ConfigurationItem` objects can
be changed at runtime by users. So you always read their values instead
of just storing their initial value to some other variable (or used as a
default for a function). For example, the following will work, but is
incorrect usage::

    SOME_SETTING = ConfigurationItem('SOME_SETTING',1,'A description.')

    def some_func(val=SOME_SETTING()):
        return val + 2

This works fine as long as the user doesn't change its value during
runtime, but if they do, the function won't know about the change::

    >>> some_func()
    3
    >>> SOME_SETTING.set(3)
    >>> some_func()  # naively should return 5, because 3 + 2 = 5
    3

There are two ways around this.  The typical/intended way is::

    def some_func():
        """
        The `SOME_SETTING` configuration item influences this output
        """
        return SOME_SETTING() + 2

Or, if the option needs to be available as a function parameter::

    def some_func(val=None):
        """
        If not specified, `val` is set by the `SOME_SETTING` configuration item.
        """
        return (SOME_SETTING() if val is None else val) + 2






See Also
========

:doc:`/logging` (overview of `astropy.logger`)


Reference/API
=============

.. automodapi:: astropy.config
    :no-inheritance-diagram:
