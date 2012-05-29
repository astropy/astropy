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

The Astropy configuration system is designed to make it easy to see what general
options are to be used

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
    populate the configuration files. Hopefully it will be true soon,
    though. The
    :func:`astropy.config.configuration._generate_all_config_items`
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

These files are only read when the relevant package is imported, though,
so any changes you make after starting an astropy session will not be
detected. You can, however, change configuration items directly. This can
be accomplished using the `set` method of a `ConfigurationItem` object as the
following example shows::

    # astropy.config.io.fits has a configuration item that sets whether or not
    # FITS extension names are case-sensitive
    >>> from astropy.io import fits
    >>> fits.EXTENSION_NAME_CASE_SENSITIVE  # Assuming you're on the default
    False
    >>> f = fits.open('somefile.fits')  # it's first extension is named "AnExt"
    >>> f[1].name
    'ANEXT'
    >>> fits.EXTENSION_NAME_CASE_SENSITIVE.set(True)
    >>> f = fits.open('somefile.fits')
    >>> f[1].name
    'AnExt'

If you want this setting to persist across later sessions, you can save this
setting by following the above code with::

    >>> fits.EXTENSION_NAME_CASE_SENSITIV.save()

Alternatively, you can simply modify the value of the
"extension_name_case_sensitive" entry in the "[io.fits]" section of
``astropy.cfg`` to True, but this would require a restart of your python
session.

For Developers
--------------
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

    SOME_OPTION = ConfigurationItem('some_option',1,'A description.')
    ANOTHER_OPTION = ConfigurationItem('annother_opt','a string val',
                            'A longer description of what this does.')

    ... implementation ...
    def some_func():
        #to get the value of these options, I might do:
        something = SOME_OPTION() + 2
        return ANOTHER_OPTION() + ' Also, I added text.'

It is highly recommended that any configuration items be placed at the
top of a module like this, as they can then be easily found when viewing
the source code and the automated tools to generate the default
configuration files can also locate these items.

There are a couple important gotchas to remember about using configuration
items in your code. First, it is tempting to do something like::

    SOME_OPTION = ConfigurationItem('some_option',1,'A description.')

    def some_func():
        return SOME_OPTION + 2  # WRONG, you wanted SOME_OPTION() + 2

but this is incorrect, because ``SOME_OPTION`` instead of
``SOME_OPTION()`` will yield a `ConfigurationItem` object, instead of the
*value* of that item (an integer, in this case).

The second point to keep in mind is that `ConfigurationItem` objects can
be changed at runtime by users. So you always read their values instead
of just storing their initial value to some other variable (or used as a
default for a function). For example, the following will work, but is
incorrect usage::

    SOME_OPTION = ConfigurationItem('some_option',1,'A description.')

    def some_func(val=SOME_OPTION()):
        return val + 2

This works fine as long as the user doesn't change it's value during
runtime, but if they do, the function won't know about the change::

    >>> some_func()
    3
    >>> SOME_OPTION.set(3)
    >>> some_func()  # naively should return 5, because 3 + 2 = 5
    3

There are two ways around this.  The typical/intended way is::

    def some_func():
        """
        The `SOME_OPTION` configuration item influences this output
        """
        return SOME_OPTION() + 2

Or, if the option needs to be available as a function parameter::

    def some_func(val=None):
        """
        If not specified, `val` is set by the `SOME_OPTION` configuration item.
        """
        return (SOME_OPTION() if None else val) + 2




.. automodapi:: astropy.config
    :no-inheritance-diagram:
