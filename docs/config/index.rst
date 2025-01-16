.. _astropy_config:

***************************************
Configuration System (`astropy.config`)
***************************************

Introduction
============

The ``astropy`` configuration system is designed to give users control of
various parameters used in ``astropy`` or affiliated packages without delving
into the source code to make those changes.

.. note::
    * Before version 4.3 the configuration file was created by default
      when importing ``astropy``. Its existence was required, which is
      no longer the case.

Getting Started
===============

The ``astropy`` configuration options are most conveniently set by modifying
the configuration file. Since ``astropy`` 4.3 you need to create this file,
whereas before it was created automatically when importing ``astropy``. The
function :func:`~astropy.config.create_config_file` creates the file with all
of the default values commented out::

    >>> from astropy.config import create_config_file
    >>> create_config_file('astropy')  # doctest: +SKIP

The exact location of this file can be obtained with
:func:`~astropy.config.get_config_dir`::

    >>> from astropy.config import get_config_dir
    >>> get_config_dir()  # doctest: +SKIP

And you should see the location of your configuration directory. The default
configuration directory is ``$HOME/.astropy/config``, but this can be
customized with :ref:`environment_variables`.

.. note::
    See :ref:`astropy_config_file` for the content of this configuration file.

Once you have found the configuration file, open it with your favorite editor.
It should have all of the sections you might want, with descriptions and the
type of value that is accepted. Feel free to edit this as you wish, and
any of these changes will be reflected when you next start ``astropy``. Or call
the :func:`~astropy.config.reload_config` function if you want to see your
changes immediately in your current ``astropy`` session::

    >>> from astropy.config import reload_config
    >>> reload_config()

.. note::
    If for whatever reason your ``$HOME/.astropy`` directory is not accessible
    (i.e., you have ``astropy`` running somehow as root but you are not the root
    user), the best solution is to set :ref:`environment_variables` pointing to
    directories you control.

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
on a ``conf`` object, or using the
:meth:`~astropy.config.ConfigNamespace.set_temp` `context_manager
<https://docs.python.org/3/reference/datamodel.html#context-managers>`_, as
long as the new value complies with the `Item Types and Validation`_.

Example
^^^^^^^

..
  EXAMPLE START
  Changing the Persistent State of Configuration Values at Runtime

Suppose there is a part of your configuration file that looks like:

.. code-block:: ini

    [utils.data]

    # URL for astropy remote data site.
    dataurl = http://data.astropy.org/

    # Time to wait for remote data query (in seconds).
    remote_timeout = 10.0

If you wish to modify the ``remote_timeout`` value, but only for some small
section of your code, then :meth:`~astropy.config.ConfigNamespace.set_temp`
takes care of resetting the value you changed when you are done using it::

    >>> from astropy.utils.data import conf
    >>> conf.remote_timeout
    10.0
    >>> # Change remote_timeout, but only inside the with-statement.
    >>> with conf.set_temp('remote_timeout', 4.5):
    ...    conf.remote_timeout
    4.5
    >>> conf.remote_timeout
    10.0

You can also modify the values at runtime directly::

    >>> conf.dataurl
    'http://data.astropy.org/'
    >>> conf.dataurl = 'http://astropydata.mywebsite.com'
    >>> conf.dataurl
    'http://astropydata.mywebsite.com'
    >>> conf.remote_timeout
    10.0
    >>> conf.remote_timeout = 4.5
    >>> conf.remote_timeout
    4.5

..
  EXAMPLE END

Reloading Configuration
-----------------------

Instead of modifying the variables in Python, you can also modify the
configuration files and then use the
:meth:`~astropy.config.ConfigNamespace.reload` method.

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

This should update the variables with the values from the configuration file:

.. doctest-skip::

    >>> conf.dataurl
    'http://myotherdata.mywebsite.com/'
    >>> conf.remote_timeout
    6.3

You can reload all configuration parameters of a ``conf`` object at
once by calling :meth:`~astropy.config.ConfigNamespace.reload` with no
parameters::

    >>> conf.reload()

Or if you want to reload all the configuration items at once, not just the ones
in the module ``conf`` belongs to, use the
:func:`~astropy.config.reload_config` function::

    >>> from astropy import config
    >>> config.reload_config()

You can also reset a configuration parameter back to its default value with the
:meth:`~astropy.config.ConfigNamespace.reset` method. Note that this is the
default value defined in the Python code, which might be different from the
value in the configuration file::

    >>> conf.reset('dataurl')
    >>> conf.dataurl
    'http://data.astropy.org/'

..
  EXAMPLE END

Exploring Configuration
-----------------------

To see what configuration parameters are defined for a given ``conf``::

    >>> from astropy.utils.iers import conf
    >>> list(conf)
    ['auto_download',
     'auto_max_age',
     ...,
     'ietf_leap_second_auto_url']

You can see more detailed information about a configuration parameter by
calling the :meth:`~astropy.config.ConfigNamespace.help` method::

    >>> conf.help("auto_max_age")
    ConfigItem: auto_max_age
      cfgtype='float'
      defaultvalue=30.0
      description='Maximum age (days) of predictive data before auto-downloading. See "Auto refresh behavior" in astropy.utils.iers documentation for details. Default is 30.'
      module=astropy.utils.iers.iers
      value=30.0

You can see information about all the configuration parameters by calling
:meth:`~astropy.config.ConfigNamespace.help` without arguments::

    >>> conf.help()
    Configuration parameters for `astropy.utils.iers`.
    <BLANKLINE>
    ConfigItem: auto_download
      cfgtype='boolean'
      ...

You can also iterate through ``conf`` in a dictionary-like fashion::

    >>> for (key, cfgitem) in conf.items():
    ...     print(f'{key} default value is {cfgitem.defaultvalue}')
    auto_download default value is True
    auto_max_age default value is 30.0
      ...

Upgrading ``astropy``
---------------------

Each time you upgrade to a new major version of ``astropy``, the
configuration parameters may have changed. If you want to create a new
configuration file, you can run::

    >>> from astropy.config import create_config_file
    >>> create_config_file('astropy', overwrite=True)  # doctest: +SKIP

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
items, but should instead be :class:`~astropy.utils.state.ScienceState`
instances, so it is possible to reproduce science results without them being
affected by configuration parameters set in a particular environment.
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
        some_list = _config.ConfigItem(
            [], 'Description of some_list', cfgtype='list')
        another_list = _config.ConfigItem(
            ['value'], 'Description of another_setting', cfgtype='list')

    # Create an instance for the user
    conf = Conf()

    # implementation ...
    def some_func():
        # to get the value of some of these options, I might do:
        something = conf.some_setting + 2
        return conf.another_setting + ' Also, I added text.'

For an affiliated package called, for example, ``packagename``, the
configuration file can be generated with the
:func:`~astropy.config.create_config_file` function like this::

    >>> from astropy.config import create_config_file
    >>> create_config_file('packagename')  # doctest: +SKIP

The following content would be written to the config file template:

.. code-block:: ini

    [subpackage]
    ## Description of some_setting
    # some_setting = 1

    ## Description of another_setting
    # another_setting = foo

    ## Description of some_list
    # some_list = ,

    ## Description of another_setting
    # another_list = value,

Note that the key/value pairs are commented out. This will allow for
changing the default values in a future version of the package without
requiring the user to edit their configuration file to take advantage
of the new defaults. By convention, the descriptions of each
parameter are in comment lines starting with two hash characters
(``##``) to distinguish them from commented out key/value pairs.

Item Types and Validation
-------------------------

If not otherwise specified, a :class:`~astropy.config.ConfigItem` gets its type
from the type of the ``defaultvalue`` it is given when it is created. The item
can only ever have a value of this type, although in some cases a provided
value can be automatically converted. For example

::

    >>> conf.auto_download
    True
    >>> conf.auto_download = 0
    >>> conf.auto_download
    False

succeeds because the :class:`int` ``0`` can be safely converted to the
:class:`bool` `False`. On the other hand

::

    >>> from astropy.utils.data import conf
    >>> conf.compute_hash_block_size
    65536
    >>> conf.compute_hash_block_size = 65536.0
    Traceback (most recent call last):
    ...
    TypeError: Provided value for configuration item compute_hash_block_size
    not valid: the value "65536.0" is of the wrong type.

fails because converting a :class:`float` to an :class:`int` can lose
information.

Note that if you want the configuration item to be limited to a particular set
of options, you should pass in a :class:`list` as the ``defaultvalue`` option.
The first entry in the :class:`list` will be taken as the default, and the
:class:`list` as a whole gives all of the valid options. For example::

    an_option = ConfigItem(['a', 'b', 'c'],
                           "This option can be 'a', 'b', or 'c'")
    conf.an_option = 'b'  # succeeds
    conf.an_option = 'c'  # succeeds
    conf.an_option = 'd'  # fails!
    conf.an_option = 6    # fails!

Finally, a :class:`~astropy.config.ConfigItem` can be explicitly given a type
via the ``cfgtype`` option::

    an_int_setting = ConfigItem(
        1, 'A description.', cfgtype='integer')
    ...
    conf.an_int_setting = 3     # works fine
    conf.an_int_setting = 4.2   # fails!

If the default value's type does not match ``cfgtype``, the
:class:`~astropy.config.ConfigItem` cannot be created.

In summary, the default behavior (of automatically determining ``cfgtype``)
is usually what you want. The main exception is when you want your
configuration item to be a :class:`list`. The default behavior will treat that
as a *list of options* unless you explicitly tell it that the
:class:`~astropy.config.ConfigItem` itself is supposed to be a :class:`list`::

    # The setting must be 1, 2 or 3
    a_list_setting = ConfigItem([1, 2, 3], 'A description.')

    # The setting must be a list and is [1, 2, 3] by default
    a_list_setting = ConfigItem([1, 2, 3], 'A description.', cfgtype='list')

Details of all of the valid ``cfgtype`` items can be found in the
`validation section of the configobj manual
<https://configobj.readthedocs.io/en/latest/validate.html#the-standard-functions>`_.
Below is a list of the valid values here for quick reference:

* ``'integer'``
* ``'float'``
* ``'boolean'``
* ``'string'``
* ``'ip_addr'``
* ``'list'``
* ``'tuple'``
* ``'int_list'``
* ``'float_list'``
* ``'bool_list'``
* ``'string_list'``
* ``'ip_addr_list'``
* ``'mixed_list'``
* ``'option'``
* ``'pass'``

Usage Tips
----------

Keep in mind that :class:`~astropy.config.ConfigItem` objects can be
changed at runtime by users. So it is always recommended to read their
values immediately before use instead of storing their initial
value to some other variable (or used as a default for a
function). For example, the following is incorrect usage::

    >>> from astropy.utils.data import conf
    >>> conf.remote_timeout = 1.0
    >>> def some_func(val=conf.remote_timeout):
    ...     return val + 2

This works only as long as the user does not change the value of the
configuration item after the function has been defined, but if they do, the
function will not know about the change::

    >>> some_func()
    3.0
    >>> conf.remote_timeout = 3.0
    >>> some_func()  # naively should return 5.0, because 3 + 2 = 5
    3.0

There are two ways around this. The typical/intended way is::

    >>> def some_func():
    ...     """
    ...     The `SOME_SETTING` configuration item influences this output
    ...     """
    ...     return conf.remote_timeout + 2
    >>> some_func()
    5.0
    >>> conf.remote_timeout = 5.0
    >>> some_func()
    7.0

Or, if the option needs to be available as a function parameter::

    def some_func(val=None):
        """
        If not specified, `val` is set by the `SOME_SETTING` configuration item.
        """
        return (conf.remote_timeout if val is None else val) + 2


Customizing Config Location in Affiliated Packages
==================================================

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

:doc:`/logging` (overview of `astropy.logger`)


Reference/API
=============

.. toctree::
   :maxdepth: 2

   ref_api

.. testcleanup::

    >>> from astropy import config
    >>> config.reload_config()
    >>> from astropy.utils.iers import conf
    >>> conf.auto_download = False
