.. currentmodule:: astropy

.. _environment_variables:

*********************
Environment variables
*********************

Cache and configuration locations
=================================

Since version 8.0.0, astropy follows the XDG specification by default, on every platform
(including non-Linux ones). As a result, the default cache location is
``$XDG_CACHE_HOME/astropy``, where ``XDG_CACHE_HOME`` itself defaults to
``$HOME/.cache``. The same goes for configuration, replacing ``cache`` with ``config``,
and preserving case.

In addition to these, and since v8.0.0, astropy supports comparable, tool-specific
environment variables for finer control:

.. glossary::

    ``ASTROPY_CACHE_DIR``
        takes precedence over ``XDG_CACHE_HOME`` and defines an
        entire path (as opposed to ``XDG_CACHE_HOME`` which only defines the *parent*
        directory of the one used by astropy). Its value must represent an absolute path,
        and must not point to file. Invalid values are ignored with a warning when the
        variable is evaluated.
        See :ref:`utils-data` for how to programmatically set or get the location of the
        corresponding directory at runtime.

    ``ASTROPY_CONFIG_DIR``
        takes precedence over ``XDG_CONFIG_HOME`` and defines an
        entire path (as opposed to ``XDG_CONFIG_HOME`` which only defines the *parent*
        directory of the one used by astropy). Its value must represent an absolute path,
        and must not point to file. Invalid values are ignored with a warning when the
        variable is evaluated.
        See :ref:`astropy_config` for how to programmatically set or get the location of
        the corresponding directory at runtime.
