.. doctest-skip-all

.. _astropy_config:

***************************************
Configuration System (`astropy.config`)
***************************************

Introduction
============

The Astropy configuration system is designed to give users control of various
parameters used in Astropy or affiliated packages without delving into the
source code to make those changes.

.. note::
    * Before version 4.3 the configuration file was created by default
      when importing ``astropy``. Its existence was required, which is
      no longer the case.
    * The configuration system got a major overhaul in ``astropy`` 0.4 as
      part of APE3. See :ref:`config-0-4-transition` for information
      about updating code to use the new API.


Getting Started
===============

The Astropy configuration options are most conveniently set by modifying the
configuration file. Since Astropy 4.3 you need to create this file, whereas
before it was created automatically when importing Astropy. To create the file
with all of the default values commented out::

    >>> from astropy.config import create_config_file
    >>> create_config_file('astropy')

You can find the exact location of this file by doing::

    >>> from astropy.config import get_config_dir
    >>> get_config_dir()

And you should see the location of your configuration directory. The standard
scheme generally puts your configuration directory in
``$HOME/.astropy/config``. It can be customized with the environment variable
``XDG_CONFIG_HOME`` and the ``$XDG_CONFIG_HOME/astropy`` directory must exist.
Note that ``XDG_CONFIG_HOME`` comes from a Linux-centric specification (see
`here <https://wiki.archlinux.org/index.php/XDG_Base_Directory_support>`_ for
more details), but Astropy will use this on any OS as a more general means to
know where user-specific configurations should be written.

.. note::
    See :ref:`astropy_config_file` for the content of this configuration file.

Once you have found the configuration file, open it with your favorite editor.
It should have all of the sections you might want, with descriptions and the
type of value that is accepted. Feel free to edit this as you wish, and
any of these changes will be reflected when you next start Astropy. Or, if you
want to see your changes immediately in your current Astropy session, do::

    >>> from astropy.config import reload_config
    >>> reload_config()

.. note::
    If for whatever reason your ``$HOME/.astropy`` directory is not accessible
    (i.e., you have ``astropy`` running somehow as root but you are not the root
    user), the best solution is to set the ``XDG_CONFIG_HOME`` and
    ``XDG_CACHE_HOME`` environment variables pointing to directories, and create
    an ``astropy`` directory inside each of those. Both the configuration and
    data download systems will then use those directories and never try to
    access the ``$HOME/.astropy`` directory.


Using `astropy.config`
======================

Accessing Values
----------------

By convention, configuration parameters live inside of objects called
``conf`` at the root of each sub-package. For example, configuration
parameters related to data files live in ``astropy.utils.data.conf``.
This object has properties for getting and setting individual
configuration parameters. For instance, to get the default URL for
``astropy`` remote data, do::

    >>> from astropy.utils.data import conf
    >>> conf.dataurl
    'http://data.astropy.org/'

Changing Values at Runtime
--------------------------

Changing the persistent state of configuration values is done by editing the
configuration file as described above. Values can also, however, be
modified in an active Python session by setting any of the properties
on a ``conf`` object.

Example
^^^^^^^

..
  EXAMPLE START
  Changing the Persistent State of Configuration Values at Runtime

If there is a part of your configuration file that looks like:

.. code-block:: ini

    [utils.data]

    # URL for astropy remote data site.
    dataurl = http://data.astropy.org/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 3.0

You should be able to modify the values at runtime this way::

    >>> from astropy.utils.data import conf
    >>> conf.dataurl
    'http://data.astropy.org/'
    >>> conf.dataurl = 'http://astropydata.mywebsite.com'
    >>> conf.dataurl
    'http://astropydata.mywebsite.com'
    >>> conf.remote_timeout
    3.0
    >>> conf.remote_timeout = 4.5
    >>> conf.remote_timeout
    4.5

..
  EXAMPLE END

Reloading Configuration
-----------------------

Instead of modifying the variables in Python, you can also modify the
configuration files and then reload them.

Example
^^^^^^^

..
  EXAMPLE START
  Modifying and Reloading Configuration Files

If you modify the configuration file to say:

.. code-block:: ini

    [utils.data]

    # URL for astropy remote data site.
    dataurl = http://myotherdata.mywebsite.com/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 6.3

And then run the following commands::

    >>> conf.reload('dataurl')
    >>> conf.reload('remote_timeout')

This should update the variables with the values from the configuration file::

    >>> conf.dataurl
    'http://myotherdata.mywebsite.com/'
    >>> conf.remote_timeout
    6.3

You can reload all configuration parameters of a ``conf`` object at
once by calling ``reload`` with no parameters::

    >>> conf.reload()

Or if you want to reload all Astropy configuration at once, use the
`~astropy.config.reload_config` function::

    >>> from astropy import config
    >>> config.reload_config('astropy')

You can also reset a configuration parameter back to its default value. Note
that this is the default value defined in the Python code, and has nothing to
do with the configuration file on disk::

    >>> conf.reset('dataurl')
    >>> conf.dataurl
    'http://data.astropy.org/'

..
  EXAMPLE END

Exploring Configuration
-----------------------

To see what configuration parameters are defined for a given ``conf``::

    >>> from astropy.utils.iers import conf
    >>> [key for key in conf]
    ['auto_download',
     'auto_max_age',
     ...,
     'ietf_leap_second_auto_url']
    >>> conf.auto_max_age
    30.0

You can also iterate through ``conf`` in a dictionary-like fashion::

    >>> [key for key in conf.keys()]
    ['auto_download',
     'auto_max_age',
     ...,
     'ietf_leap_second_auto_url']
    >>> [cfgitem for cfgitem in conf.values()]
    [<ConfigItem: name='auto_download' value=True at ...>,
     <ConfigItem: name='auto_max_age' value=30.0 at ...>,
     ...,
     <ConfigItem: name='ietf_leap_second_auto_url' value=...>]
    >>> for (key, cfgitem) in conf.items():
    ...     if key == 'auto_max_age':
    ...         print(f'{cfgitem.description} Value is {cfgitem()}')
    Maximum age (days) of predictive data before auto-downloading. Default is 30. Value is 30.0

Upgrading ``astropy``
---------------------

Each time you upgrade to a new major version of ``astropy``, the
configuration parameters may have changed. If you want to create a new
configuration file, you can run::

    >>> from astropy.config import create_config_file
    >>> create_config_file('astropy', overwrite=True)

Note that this will overwrite the existing file, so if you modified it you
may want to report your changes in the new file. Another possibility is to
have a look at the :ref:`astropy_config` to see what has changed.

.. _config-developer:

Adding New Configuration Items
==============================

Configuration items should be used wherever an option or setting is
needed that is either tied to a system configuration or should persist
across sessions of ``astropy`` or an affiliated package. Options that may
affect the results of science calculations should not be configuration
items, but should instead be `astropy.utils.state.ScienceState`, so
it is possible to reproduce science results without them being affected
by configuration parameters set in a particular environment.
Admittedly, this is only a guideline, as the precise cases where a
configuration item is preferred over, say, a keyword option for a
function is somewhat personal preference. It is the preferred form of
persistent configuration, however, and ``astropy`` packages must all use
it (and it is recommended for affiliated packages).

The reference guide below describes the interface for creating a
``conf`` object with a number of configuration parameters. They
should be defined at the top level (i.e., in the ``__init__.py`` of
each sub-package that has configuration items)::

    """ This is the docstring at the beginning of a module
    """
    from astropy import config as _config

    class Conf(_config.ConfigNamespace):
        """
        Configuration parameters for my subpackage.
        """
        some_setting = _config.ConfigItem(
            1, 'Description of some_setting')
        another_setting = _config.ConfigItem(
            'string value', 'Description of another_setting')
    # Create an instance for the user
    conf = Conf()

    ... implementation ...
    def some_func():
        #to get the value of these options, I might do:
        something = conf.some_setting + 2
        return conf.another_setting + ' Also, I added text.'

The configuration items also need to be added to the config file
template. For ``astropy``, this file is in ``astropy/astropy.cfg``. For
an affiliated package called, for example, ``packagename``, the file
is in ``packagename/packagename.cfg``. For the example above, the
following content would be added to the config file template:

.. code-block:: ini

    [subpackage]
    ## Description of some_setting
    # some_setting = 1

    ## Description of another_setting
    # another_setting = foo

Note that the key/value pairs are commented out. This will allow for
changing the default values in a future version of ``astropy`` without
requiring the user to edit their configuration file to take advantage
of the new defaults. By convention, the descriptions of each
parameter are in comment lines starting with two hash characters
(``##``) to distinguish them from commented out key/value pairs.

Item Types and Validation
-------------------------

If not otherwise specified, a `~astropy.config.ConfigItem` gets its
type from the type of the ``defaultvalue`` it is given when it is
created. The item can only be set to be an object of this type.
Hence::

    some_setting = ConfigItem(1, 'A description.')
    ...
    conf.some_setting = 1.2

will fail, because ``1.2`` is a float and ``1`` is an integer.

Note that if you want the configuration item to be limited to a
particular set of options, you should pass in a list as the
``defaultvalue`` option. The first entry in the list will be taken as
the default, and the list as a whole gives all of the valid options. For
example::

    an_option = ConfigItem(
        ['a', 'b', 'c'],
        "This option can be 'a', 'b', or 'c'")
    ...
    conf.an_option = 'b'  # succeeds
    conf.an_option = 'c'  # succeeds
    conf.an_option = 'd'  # fails!
    conf.an_option = 6    # fails!

Finally, a `~astropy.config.ConfigItem` can be explicitly given a type
via the ``cfgtype`` option::

    an_int_setting = ConfigItem(
        1, 'A description.', cfgtype='integer')
    ...
    conf.an_int_setting = 3     # works fine
    conf.an_int_setting = 4.2   # fails!

If the default value's type does not match ``cfgtype``, the
`~astropy.config.ConfigItem` cannot be created::

    an_int_setting = ConfigItem(
        4.2, 'A description.', cfgtype='integer')

In summary, the default behavior (of automatically determining ``cfgtype``)
is usually what you want. The main exception is when you want your
configuration item to be a list. The default behavior will treat that
as a list of *options* unless you explicitly tell it that the
`~astropy.config.ConfigItem` itself is supposed to be a list::

    a_list_setting = ConfigItem([1, 2, 3], 'A description.')

    a_list_setting = ConfigItem([1, 2, 3], 'A description.', cfgtype='list')

Details of all of the valid ``cfgtype`` items can be found in the
`validation section of the configobj manual
<http://www.voidspace.org.uk/python/validate.html#the-standard-functions>`_.
Below is a list of the valid values here for quick reference:

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
----------

Keep in mind that `~astropy.config.ConfigItem` objects can be
changed at runtime by users. So it is always recommended to read their
values immediately before use instead of storing their initial
value to some other variable (or used as a default for a
function). For example, the following will work, but is incorrect
usage::

    def some_func(val=conf.some_setting):
        return val + 2

This works fine as long as the user does not change its value during
runtime, but if they do, the function will not know about the change::

    >>> some_func()
    3
    >>> conf.some_setting = 3
    >>> some_func()  # naively should return 5, because 3 + 2 = 5
    3

There are two ways around this. The typical/intended way is::

    def some_func():
        """
        The `SOME_SETTING` configuration item influences this output
        """
        return conf.some_setting + 2

Or, if the option needs to be available as a function parameter::

    def some_func(val=None):
        """
        If not specified, `val` is set by the `SOME_SETTING` configuration item.
        """
        return (conf.some_setting if val is None else val) + 2


Customizing Config in Affiliated Packages
=========================================

The `astropy.config` package can be used by other packages. By default creating
a config object in another package will lead to a configuration file taking the
name of that package in the ``astropy`` config directory (i.e.,
``<astropy_config>/packagename.cfg``).


It is possible to configure this behavior so that the a custom configuration
directory is created for your package, for example,
``~/.packagename/packagename.cfg``. To do this, create a ``packagename.config``
subpackage and put the following into the ``__init__.py`` file::

  import astropy.config as astropyconfig


  class ConfigNamespace(astropyconfig.ConfigNamespace):
      rootname = 'packagename'


  class ConfigItem(astropyconfig.ConfigItem):
      rootname = 'packagename'

Then replace all imports of `astropy.config` with ``packagename.config``.


See Also
========

.. toctree::
   :maxdepth: 2

   astropy_config
   config_0_4_transition

:doc:`/logging` (overview of `astropy.logger`)


Reference/API
=============

.. automodapi:: astropy.config
