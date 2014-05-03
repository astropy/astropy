.. doctest-skip-all

.. _astropy_config:

***************************************
Configuration system (`astropy.config`)
***************************************

Introduction
============

The astropy configuration system is designed to give users control of various
parameters used in astropy or affiliated packages without delving into the
source code to make those changes.

.. note::
    The configuration system got a major overhaul in astropy 0.4 as
    part of APE3.  See :ref:`config-0-4-transition` for information
    about updating code to use the new API.


Getting Started
===============

The Astropy configuration options are most easily set by modifying the
configuration file.  It will be automatically generated with all the
default values commented out the first time you import Astropy.  You
can find the exact location by doing::

    >>> from astropy.config import get_config_dir
    >>> get_config_dir()

And you should see the location of your configuration directory. The standard
scheme generally puts your configuration directory in ``$HOME/.astropy/config``,
but if you've set the environment variable ``XDG_CONFIG_HOME`` and the
``$XDG_CONFIG_HOME/astropy`` directory exists, it will instead be there.

.. note::
    This is a slight variation from the behavior of most Linux
    applications with respect to ``$XDG_CONFIG_HOME``, where the
    default, if ``$XDG_CONFIG_HOME`` is not defined, would be to put
    configuration in ``~/.config/astropy``.

Once you've found the configuration file, open it with your favorite editor.
It should have all of the sections you might want, with descriptions and the
type of the value that is accepted.  Feel free to edit this as you wish, and
any of these changes will be reflected when you next start Astropy.  Or, if you
want to see your changes immediately in your current Astropy session, just do::

    >>> from astropy.config import reload_config
    >>> reload_config()

.. note::
    If for whatever reason your ``$HOME/.astropy`` directory is not accessible
    (i.e., you have astropy running somehow as root but you are not the root
    user), the best solution is to set the ``XDG_CONFIG_HOME`` and
    ``XDG_CACHE_HOME`` environment variables pointing to directories, and create
    an ``astropy`` directory inside each of those.  Both the configuration and
    data download systems will then use those directories and never try to
    access the ``$HOME/.astropy`` directory.


Using `astropy.config`
======================

Accessing Values
----------------

By convention, configuration parameters live inside of objects called
``conf`` at the root of each subpackage.  For example, configuration
parameters related to data files live in ``astropy.utils.data.conf``.
This object has properties for getting and setting individual
configuration parameters.  For instance to get the default URL for
astropy remote data do::

    >>> from astropy.utils.data import conf
    >>> conf.dataurl
    'http://data.astropy.org/'

Changing Values at Run-time
---------------------------

Changing configuration values persistently is done by editing the
configuration file as described above.  Values can also, however, be
modified in an active python session by setting any of the properties
on a ``conf`` object.

For example, if there is a part of your configuration file that looks
like:

.. code-block:: ini

    [utils.data]

    # URL for astropy remote data site.
    dataurl = http://data.astropy.org/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 3.0

You should be able to modify the values at run-time this way::

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

Reloading Configuration
-----------------------

Instead of modifying the variables in python, you can also modify the
configuration files and then reload them.  For example, if you modify the
configuration file to say:

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

Or if you want to reload all astropy configuration at once, use the
`~astropy.config.reload_config` function::

    >>> config.reload_config('astropy')

You can also reset a configuration parameter back to its default value.  Note that this is the default value defined in the Python code, and has nothing to do with the configuration file on disk::

    >>> conf.reset('dataurl')
    >>> conf.dataurl
    'http://data.astropy.org/'

Upgrading astropy
-----------------

Each time you upgrade to a new major version of astropy, the
configuration parameters may have changed.

If you never edited your configuration file, there is nothing for you
to do.  It will automatically be replaced with a configuration file
template for the newly installed version of astropy.

If you did customize your configuration file, it will not be touched.
Instead, a new configuration file template will be installed alongside
it with the version number in the filename, for example
``astropy.0.4.cfg``.  You can compare this file to your
``astropy.cfg`` file to see what needs to be changed or updated.

.. _config-developer:

Adding new configuration items
==============================

Configuration items should be used wherever an option or setting is
needed that is either tied to a system configuration or should persist
across sessions of astropy or an affiliated package.  Options that may
affect the results of science calculations should not be configuration
items, but should instead be `astropy.utils.state.ScienceState`, so
it's possible to reproduce science results without them being affected
by configuration parameters set in a particular environment.
Admittedly, this is only a guideline, as the precise cases where a
configuration item is preferred over, say, a keyword option for a
function is somewhat personal preference. It is the preferred form of
persistent configuration, however, and astropy packages must all use
it (and it is recommended for affiliated packages).

The reference guide below describes the interface for creating a
``conf`` object with a number of configuration parameters.  They
should be defined at the top level, i.e. in the ``__init__.py`` of
each subpackage that has configuration items::

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
template.  For astropy, this file is in ``astropy/astropy.cfg``.  For
an affiliated package called, for example, ``packagename``, the file
is in ``packagename/packagename.cfg``.  For the example above, the
following content would be added to the config file template:

.. code-block:: ini

    [subpackage]
    ## Description of some_setting
    # some_setting = 1

    ## Description of another_setting
    # another_setting = foo

Note that the key/value pairs are commented out.  This will allow for
changing the default values in a future version of astropy without
requiring the user to edit their configuration file to take advantage
of the new defaults.  By convention, the descriptions of each
parameter are in comment lines starting with two hash characters
(``##``) to distinguish them from commented out key/value pairs.

Item Types and Validation
-------------------------

If not otherwise specified, a `~astropy.config.ConfigItem` gets its
type from the type of the ``defaultvalue`` it is given when it is
created.  The item can only be set to be an object of this type.
Hence::

    some_setting = ConfigItem(1, 'A description.')
    ...
    conf.some_setting = 1.2

will fail, because ``1.2`` is a float and ``1`` is an int.

Note that if you want the configuration item to be limited to a
particular set of options, you should pass in a list as the
``defaultvalue`` option.  The first entry in the list will be taken as
the default, and the list as a whole gives all the valid options.  For
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

If the default value's type doesn't match ``cfgtype``, the
`~astropy.config.ConfigItem` cannot be created::

    an_int_setting = ConfigItem(
        4.2, 'A description.', cfgtype='integer')

In summary, the default behavior (of automatically determining ``cfgtype``)
is usually what you want.  The main exception is when you want your
configuration item to be a list.  The default behavior will treat that
as a list of *options* unless you explicitly tell it that the
`~astropy.config.ConfigItem` itself is supposed to be a list::

    a_list_setting = ConfigItem([1, 2, 3], 'A description.')

    a_list_setting = ConfigItem([1, 2, 3], 'A description.', cfgtype='list')

Details of all the valid ``cfgtype`` items can be found in the
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

Keep in mind is that `~astropy.config.ConfigItem` objects can be
changed at runtime by users. So it is always recommended to read their
values immediately before use instead of just storing their initial
value to some other variable (or used as a default for a
function). For example, the following will work, but is incorrect
usage::

    def some_func(val=conf.some_setting):
        return val + 2

This works fine as long as the user doesn't change its value during
runtime, but if they do, the function won't know about the change::

    >>> some_func()
    3
    >>> conf.some_setting = 3
    >>> some_func()  # naively should return 5, because 3 + 2 = 5
    3

There are two ways around this.  The typical/intended way is::

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


See Also
========

.. toctree::
   :maxdepth: 2

   config_0_4_transition

:doc:`/logging` (overview of `astropy.logger`)


Reference/API
=============

.. automodapi:: astropy.config
