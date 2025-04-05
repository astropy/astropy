.. currentmodule:: astropy

.. _environment_variables:

*********************
Environment variables
*********************


Runtime
-------

.. glossary::

    ``ASTROPY_CONFIG_DIR`` / ``XDG_CONFIG_HOME``
        These environment variables control where configuration files are read
        and written. Astropy will look for configuration files in
        ``$ASTROPY_CONFIG_DIR/astropy``, respectively
        ``$XDG_CONFIG_HOME/astropy``, depending on which variable is defined,
        the former taking priority.
        See :ref:`astropy_config` for how to programmatically set or get the
        location of the corresponding directory at runtime.
        ``ASTROPY_CONFIG_DIR`` was introduced in astropy 7.1 and as such, it
        will be ignored by older versions.

    ``ASTROPY_CACHE_DIR`` / ``XDG_CACHE_HOME``
        These environment variables control where data files are cached.
        Astropy will cache files in ``$ASTROPY_CACHE_DIR/astropy``, respectively
        ``$XDG_CACHE_HOME/astropy``, depending on which variable is defined,
        the former taking priority.
        See :ref:`utils-data` for how to programmatically set or get the
        location of the corresponding directory at runtime.
        ``ASTROPY_CACHE_DIR`` was introduced in astropy 7.1 and as such, it
        will be ignored by older versions.


.. note::
    ``XDG_CONFIG_HOME`` and ``XDG_CACHE_HOME`` come from a Linux-centric
    specification
    (see `here <https://wiki.archlinux.org/index.php/XDG_Base_Directory_support>`_
    for more details), but ``astropy`` will use them on any OS as a more
    general means to know where user-specific configurations should be written.
