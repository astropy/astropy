# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Utility functions and classes"""

from __future__ import print_function, division

import inspect
import traceback

from ...extern.six.moves import queue
from ...extern.six.moves.urllib.error import URLError
from ...extern.six.moves.urllib.request import urlopen
from ...extern.six import StringIO
from ...config import ConfigurationItem

from .constants import SAMP_STATUS_ERROR

ALLOW_INTERNET = ConfigurationItem('use_internet', True,
                                   "Whether to allow astropy.vo.samp to use the internet, if available")


def internet_on():
    if not ALLOW_INTERNET():
        return False
    else:
        try:
            urlopen('http://google.com', timeout=1)
            return True
        except URLError:
            pass
        return False

__all__ = ["SAMPMsgReplierWrapper"]

__doctest_skip__ = ['.']


class _ServerProxyPoolMethod:
    # some magic to bind an XML-RPC method to an RPC server.
    # supports "nested" methods (e.g. examples.getStateName)

    def __init__(self, proxies, name):
        self.__proxies = proxies
        self.__name = name

    def __getattr__(self, name):
        return _ServerProxyPoolMethod(self.__proxies, "%s.%s" % (self.__name, name))

    def __call__(self, *args, **kwrds):
        proxy = self.__proxies.get()
        try:
            response = eval("proxy.%s(*args, **kwrds)" % self.__name)
        except:
            self.__proxies.put(proxy)
            raise
        self.__proxies.put(proxy)
        return response


class ServerProxyPool(object):
    """
    A thread-safe pool of `xmlrpc.ServerProxy` objects.
    """

    def __init__(self, size, proxy_class, *args, **keywords):

        self._proxies = queue.Queue(size)
        for i in range(size):
            self._proxies.put(proxy_class(*args, **keywords))

    def __getattr__(self, name):
        # magic method dispatcher
        return _ServerProxyPoolMethod(self._proxies, name)


class SAMPMsgReplierWrapper(object):
    """
    Decorator class/function that allows to automatically grab
    errors and returned maps (if any) from a function bound
    to a SAMP call (or notify).

    Parameters
    ----------
    cli : `SAMPIntegratedClient` or `SAMPClient`
        SAMP client instance.
        Decorator initialization, accepting the instance of the
        client that receives the call or notification.
    """

    def __init__(self, cli):
        self.cli = cli

    def __call__(self, f):

        def wrapped_f(*args):

            if ((inspect.ismethod(f) and f.__func__.__code__.co_argcount == 6)
                or (inspect.isfunction(f) and f.__code__.co_argcount == 5)
                    or args[2] is None):

                # It is a notification
                f(*args)

            else:
                # It's a call
                try:
                    result = f(*args)
                    if result:
                        self.cli.hub.reply(self.cli.get_private_key(), args[2],
                                           {"samp.status": SAMP_STATUS_ERROR,
                                            "samp.result": result})
                except:
                    err = StringIO()
                    traceback.print_exc(file=err)
                    txt = err.getvalue()
                    self.cli.hub.reply(self.cli.get_private_key(), args[2],
                                       {"samp.status": SAMP_STATUS_ERROR,
                                        "samp.result": {"txt": txt}})

        return wrapped_f


class _HubAsClient(object):

    def __init__(self, handler):
        self._handler = handler

    def __getattr__(self, name):
        # magic method dispatcher
        return _HubAsClientMethod(self._handler, name)


class _HubAsClientMethod(object):

    def __init__(self, send, name):
        self.__send = send
        self.__name = name

    def __getattr__(self, name):
        return _HubAsClientMethod(self.__send, "%s.%s" % (self.__name, name))

    def __call__(self, *args):
        return self.__send(self.__name, args)
