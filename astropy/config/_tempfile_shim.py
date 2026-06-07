__all__ = ["TemporaryDirectory"]

import sys

if sys.version_info >= (3, 12, 1):
    from tempfile import TemporaryDirectory
else:
    import os as _os
    import shutil as _shutil
    import types as _types
    import warnings as _warnings
    import weakref as _weakref
    from tempfile import mkdtemp as _mkdtemp

    # backport the class from Python 3.12.1 for the delete argument
    # the 3.12.0 version had a couple defects, but there wasn't
    # any changes between 3.12.1 and 3.12.13, and no bugfixes between 3.12.1 and 3.14.4
    # this version should be perfectly stable on 3.11.0 to 3.12.0 (included)

    def _os_path_isjunction(path):
        """Test whether a path is a junction"""
        # Junctions are not a part of posix semantics
        is_posix = "posix" in sys.builtin_module_names
        if not is_posix and hasattr(_os.stat_result, "st_reparse_tag"):
            import stat

            try:
                st = _os.lstat(path)
            except (OSError, ValueError, AttributeError):
                return False
            return bool(st.st_reparse_tag == stat.IO_REPARSE_TAG_MOUNT_POINT)
        else:
            _os.fspath(path)
            return False

    def _dont_follow_symlinks(func, path, *args):
        # Pass follow_symlinks=False, unless not supported on this platform.
        if func in _os.supports_follow_symlinks:
            func(path, *args, follow_symlinks=False)
        elif _os.name == "nt" or not _os.path.islink(path):
            func(path, *args)

    def _resetperms(path):
        try:
            chflags = _os.chflags
        except AttributeError:
            pass
        else:
            _dont_follow_symlinks(chflags, path, 0)
        _dont_follow_symlinks(_os.chmod, path, 0o700)

    class TemporaryDirectory:
        """Create and return a temporary directory.  This has the same
        behavior as mkdtemp but can be used as a context manager.  For
        example:

            with TemporaryDirectory() as tmpdir:
                ...

        Upon exiting the context, the directory and everything contained
        in it are removed (unless delete=False is passed or an exception
        is raised during cleanup and ignore_cleanup_errors is not True).

        Optional Arguments:
            suffix - A str suffix for the directory name.  (see mkdtemp)
            prefix - A str prefix for the directory name.  (see mkdtemp)
            dir - A directory to create this temp dir in.  (see mkdtemp)
            ignore_cleanup_errors - False; ignore exceptions during cleanup?
            delete - True; whether the directory is automatically deleted.
        """

        def __init__(
            self,
            suffix=None,
            prefix=None,
            dir=None,
            ignore_cleanup_errors=False,
            *,
            delete=True,
        ):
            self.name = _mkdtemp(suffix, prefix, dir)
            self._ignore_cleanup_errors = ignore_cleanup_errors
            self._delete = delete
            self._finalizer = _weakref.finalize(
                self,
                self._cleanup,
                self.name,
                warn_message="Implicitly cleaning up {!r}".format(self),
                ignore_errors=self._ignore_cleanup_errors,
                delete=self._delete,
            )

        @classmethod
        def _rmtree(cls, name, ignore_errors=False, repeated=False):
            def onexc(func, path, exc):
                if isinstance(exc, PermissionError):
                    if repeated and path == name:
                        if ignore_errors:
                            return
                        raise

                    try:
                        if path != name:
                            _resetperms(_os.path.dirname(path))
                        _resetperms(path)

                        try:
                            _os.unlink(path)
                        except IsADirectoryError:
                            cls._rmtree(path, ignore_errors=ignore_errors)
                        except PermissionError:
                            # The PermissionError handler was originally added for
                            # FreeBSD in directories, but it seems that it is raised
                            # on Windows too.
                            # bpo-43153: Calling _rmtree again may
                            # raise NotADirectoryError and mask the PermissionError.
                            # So we must re-raise the current PermissionError if
                            # path is not a directory.
                            if not _os.path.isdir(path) or _os_path_isjunction(path):
                                if ignore_errors:
                                    return
                                raise
                            cls._rmtree(
                                path,
                                ignore_errors=ignore_errors,
                                repeated=(path == name),
                            )
                    except FileNotFoundError:
                        pass
                elif isinstance(exc, FileNotFoundError):
                    pass
                else:
                    if not ignore_errors:
                        raise

            if sys.version_info < (3, 12):
                onexc_kwarg = "onerror"
            else:
                onexc_kwarg = "onexc"

            _shutil.rmtree(name, **{onexc_kwarg: onexc})

        @classmethod
        def _cleanup(cls, name, warn_message, ignore_errors=False, delete=True):
            if delete:
                cls._rmtree(name, ignore_errors=ignore_errors)
                _warnings.warn(warn_message, ResourceWarning)

        def __repr__(self):
            return "<{} {!r}>".format(self.__class__.__name__, self.name)

        def __enter__(self):
            return self.name

        def __exit__(self, exc, value, tb):
            if self._delete:
                self.cleanup()

        def cleanup(self):
            if self._finalizer.detach() or _os.path.exists(self.name):
                self._rmtree(self.name, ignore_errors=self._ignore_cleanup_errors)

        __class_getitem__ = classmethod(_types.GenericAlias)
