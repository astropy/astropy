"""
A web cache class with the following features:

   - The network is only hit when the cached file is more than X old.

   - After every network check, the timestamp on the cached file is
     updated, so the next network check will be at least after X time.

   - N (default 5) backup copies of the cached files are maintained,
     so if the network copies are damaged, local users can restore
     back to a working state.
"""

from __future__ import division, absolute_import

# STDLIB
import datetime
import email.utils
import os
import re
import shutil
import sys
import time
import urllib2
import warnings

# LOCAL
from . import voexceptions
from .voexceptions import vo_warn, W23
from . import webquery


def parse_http_date(date):
    pass


class Cache:
    def __init__(self, local_root, remote_root, cache_delta=None,
                 keep_backup_copies=5):
        """
        A cache repository to request static resources from a web
        directory and a locally-cached mirror directory.

        Parameters
        ----------

        local_root : string path
            Location of the cache directory on the local filesystem.

        remote_root : string URL
            URL to the root of the remote databases.

        cache_delta : datetime.timedelta
            Minimum age of the local file before making a remote query
            to attempt to update the database.

        keep_backup_copies : int
            Number of backup copies of each local file to keep.
        """
        self._local_root = local_root
        self._remote_root = remote_root
        if cache_delta is None:
            cache_delta = datetime.timedelta(days=1)
        else:
            assert isinstance(cache_delta, datetime.timedelta)
        self._cache_delta = cache_delta

        cache_dir = self.get_cache_dir()
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        # TODO: Ensure writeable

        self._keep_backup_copies = int(keep_backup_copies)

        self._lock = 0

    def get_local_root_dir(self):
        """
        Get the local root directory.  If the string contains `~` it
        will be expanded to the user's home directory.
        """
        return os.path.expanduser(self._local_root)

    def get_cache_dir(self):
        """
        Get the cache directory, within the local root directory.
        """
        return os.path.join(self.get_local_root_dir(), 'cache')

    def _get_lock_dir(self):
        """
        Get the path to the lock directory
        """
        return os.path.join(self.get_local_root_dir(), 'lock')

    def _acquire_lock(self):
        """
        Try to acquire a lock on the cache directory.  This should be
        done before any action that will update/modify the cache.  May
        be called recursively.

        (This is a single lock for the whole cache -- in the future,
        it may be sufficient to lock on a per-file basis.)
        """
        self._lock += 1
        if self._lock == 1:
            lock_dir = self._get_lock_dir()
            for i in range(5):
                try:
                    os.mkdir(lock_dir)
                except OSError:
                    time.sleep(1)
                else:
                    return
            raise RuntimeError(
                "Unable to acquire lock for cache directory ('%s' exists)" %
                lock_dir)

    def _release_lock(self):
        """
        Release the lock on the cache directory.  This should always
        be called in pairs with `acquire_lock`.
        """
        self._lock -= 1
        if self._lock == 0:
            lock_dir = self._get_lock_dir()
            if os.path.exists(lock_dir) and os.path.isdir(lock_dir):
                os.removedirs(lock_dir)
            else:
                raise RuntimeError(
                    ("Error releasing lock. '%s' either does not exist or is " +
                     "not a directory.") %
                    lock_dir)

    def get_cache_file_path(self, filename):
        """
        Get the path to a file from the local cache.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.

        Returns
        -------
        path : str
            The full path to the file.
        """
        cache_dir = self.get_cache_dir()
        cache_file = os.path.abspath(os.path.join(cache_dir, filename))
        if not cache_file.startswith(cache_dir):
            raise ValueError(
                "Requested cached file outside of cache directory tree.")
        return cache_file

    def get_cache_date(self, filename):
        """
        Return the date of last modification of the file in the cache.
        Returns `0` if the file doesn't exist.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.

        Returns
        -------
        timestamp : datetime.datetime object
        """
        cache_file = self.get_cache_file_path(filename)
        if not os.path.exists(cache_file):
            return datetime.datetime.fromtimestamp(0)
        return datetime.datetime.fromtimestamp(os.stat(cache_file).st_mtime)

    def has_cached_copy(self, filename):
        """
        Returns `True` if a cached copy of the resource exists.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.
        """
        cache_file = self.get_cache_file_path(filename)
        return os.path.exists(cache_file) and os.path.isfile(cache_file)

    def get_cache_content(self, filename):
        """
        Returns a read-only file-like object to read content from the
        file.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.

        Returns
        -------
        input_file : read-only file-like object
        """
        cache_file = self.get_cache_file_path(filename)
        return open(cache_file, 'rb')

    def get_cache_backups(self, filename):
        """
        Returns a list of the backup cached files for the given file.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.

        Returns
        -------
        paths : list of paths
        """
        cache_file = self.get_cache_file_path(filename)
        dirname, fname = os.path.split(cache_file)
        candidate_files = []
        regex = re.compile(
            u"^" + fname + ur"~\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d$")
        for existing_file in os.listdir(dirname):
            if regex.match(existing_file):
                existing_file_path = os.path.join(dirname, existing_file)
                candidate_files.append((
                        datetime.datetime.fromtimestamp(
                            os.stat(existing_file_path).st_mtime),
                        existing_file_path))
        candidate_files.sort()

        return candidate_files

    def backup_cache_file(self, filename):
        """
        Backup the cached file, leaving at most
        *self.keep_backup_copies* around on disk.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.
        """
        cache_file = self.get_cache_file_path(filename)

        self._acquire_lock()
        try:
            candidate_files = self.get_cache_backups(filename)

            # Remove the oldest, keeping self.keep_backup_copies
            for date, existing_file_path in candidate_files[:-self._keep_backup_copies]:
                os.remove(os.path.join(existing_file_path))

            if os.path.exists(cache_file):
                # Create the new backup
                date = datetime.datetime.fromtimestamp(
                    os.stat(cache_file).st_ctime)
                backup_file = cache_file + "~" + date.isoformat()
                os.rename(cache_file, backup_file)
        finally:
            self._release_lock()

    def put_cache_content(self, filename, in_fd):
        """
        Read from an input stream and write it into a cache file.

        If the cache file already exists, it will be backed up to a
        file with the same name + '~'.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.

        in_fd : readable file-like object
            Object to read the file content from.
        """
        self._acquire_lock()
        try:
            cache_file = self.get_cache_file_path(filename)
            cache_dir = os.path.dirname(cache_file)
            if not os.path.exists(cache_dir):
                os.makedirs(cache_dir)

            self.backup_cache_file(filename)

            out_fd = open(cache_file, 'wb')
            try:
                shutil.copyfileobj(in_fd, out_fd, 102400)
            finally:
                out_fd.close()
        finally:
            self._release_lock()

    def touch_cache_content(self, filename):
        """
        Update the timestamp on the cached content with the last time
        the remote file was compared.

        Parameters
        ----------
        filename : str
            The name of a file in the local cache.
        """
        # Not certain the lock is required here since worst case here
        # you may occasionally have to hit the webserver more than
        # necessary.
        cache_file = self.get_cache_file_path(filename)
        self._acquire_lock()
        try:
            os.utime(cache_file, None)
        finally:
            self._release_lock()

    def clear_cache(self):
        """
        Empty the entire cache.
        """
        cache_dir = self.get_cache_dir()
        self._acquire_lock()
        try:
            shutil.removedirs(cache_dir)
        finally:
            self._release_lock()

    def get_web_date(self, response):
        """
        Return the modification time from an HTTP Response as a
        datetime object.

        Parameters
        ----------
        response : HTTPResponse

        Returns
        -------
        last_modified : datetime.datetime
        """
        date_string = response.info()['Last-Modified']
        date = datetime.datetime(*email.utils.parsedate(date_string)[:6])
        return date

    def get_file(self, filename):
        """
        Get a file from the cache.  If the cache is too old, it will
        be fetched remotely and the cache updated.

        Parameters
        ----------
        filename : str filename

        Returns
        -------
        input_file : read-only file-like object
        """
        # TODO: Implement HTTP cache control per the spec.  Currently,
        # we only use the 'Last-Modified' HTTP header and have our own
        # very bandwidth conservative cache update policy.

        cache_date = self.get_cache_date(filename)

        # We only want to hit the webserver if the cached copy is more
        # than a day old.
        if datetime.datetime.now() - cache_date > self._cache_delta:
            try:
                response = webquery.webget_open(self._remote_root + filename)
                try:
                    web_date = self.get_web_date(response)
                    # If the web file is newer, copy it to our cache ...
                    if web_date > cache_date:
                        self.put_cache_content(filename, response)
                    else:
                        # ... otherwise, touch our cached file so we won't
                        # check again for another day.
                        self.touch_cache_content(filename)
                finally:
                    response.close()
            except (urllib2.HTTPError, urllib2.URLError), e:
                # The web request can fail in all kinds of ways, but
                # we don't want to propagate the error unless we don't
                # have a cached copy to return.
                if self.has_cached_copy(filename):
                    vo_warn(W23, filename)
                else:
                    # Pass the exception on up, since there is no
                    # cached copy to help us.
                    raise e

        return self.get_cache_content(filename)
