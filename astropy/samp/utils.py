# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utility functions and classes
"""


import queue
import inspect
import traceback
from io import StringIO
import xmlrpc.client as xmlrpc
from urllib.request import urlopen

from .constants import SAMP_STATUS_ERROR
from .errors import SAMPProxyError


def internet_on():
    from . import conf
    if not conf.use_internet:
        return False
    else:
        try:
            urlopen('http://google.com', timeout=1.)
            return True
        except Exception:
            return False


__all__ = ["SAMPMsgReplierWrapper"]

__doctest_skip__ = ['.']


def getattr_recursive(variable, attribute):
    """
    Get attributes recursively.
    """
    if '.' in attribute:
        top, remaining = attribute.split('.', 1)
        return getattr_recursive(getattr(variable, top), remaining)
    else:
        return getattr(variable, attribute)


class _ServerProxyPoolMethod:

    # some magic to bind an XML-RPC method to an RPC server.
    # supports "nested" methods (e.g. examples.getStateName)

    def __init__(self, proxies, name):
        self.__proxies = proxies
        self.__name = name

    def __getattr__(self, name):
        return _ServerProxyPoolMethod(self.__proxies, f"{self.__name}.{name}")

    def __call__(self, *args, **kwrds):
        proxy = self.__proxies.get()
        function = getattr_recursive(proxy, self.__name)
        try:
            response = function(*args, **kwrds)
        except xmlrpc.Fault as exc:
            raise SAMPProxyError(exc.faultCode, exc.faultString)
        finally:
            self.__proxies.put(proxy)
        return response


class ServerProxyPool:
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


class SAMPMsgReplierWrapper:
    """
    Function decorator that allows to automatically grab errors and returned
    maps (if any) from a function bound to a SAMP call (or notify).

    Parameters
    ----------
    cli : :class:`~astropy.samp.SAMPIntegratedClient` or :class:`~astropy.samp.SAMPClient`
        SAMP client instance. Decorator initialization, accepting the instance
        of the client that receives the call or notification.
    """

    def __init__(self, cli):
        self.cli = cli

    def __call__(self, f):

        def wrapped_f(*args):

            if get_num_args(f) == 5 or args[2] is None:  # notification

                f(*args)

            else:  # call

                try:
                    result = f(*args)
                    if result:
                        self.cli.hub.reply(self.cli.get_private_key(), args[2],
                                           {"samp.status": SAMP_STATUS_ERROR,
                                            "samp.result": result})
                except Exception:
                    err = StringIO()
                    traceback.print_exc(file=err)
                    txt = err.getvalue()
                    self.cli.hub.reply(self.cli.get_private_key(), args[2],
                                       {"samp.status": SAMP_STATUS_ERROR,
                                        "samp.result": {"txt": txt}})

        return wrapped_f


class _HubAsClient:

    def __init__(self, handler):
        self._handler = handler

    def __getattr__(self, name):
        # magic method dispatcher
        return _HubAsClientMethod(self._handler, name)


class _HubAsClientMethod:

    def __init__(self, send, name):
        self.__send = send
        self.__name = name

    def __getattr__(self, name):
        return _HubAsClientMethod(self.__send, f"{self.__name}.{name}")

    def __call__(self, *args):
        return self.__send(self.__name, args)


def get_num_args(f):
    """
    Find the number of arguments a function or method takes (excluding ``self``).
    """
    if inspect.ismethod(f):
        return f.__func__.__code__.co_argcount - 1
    elif inspect.isfunction(f):
        return f.__code__.co_argcount
    else:
        raise TypeError("f should be a function or a method")
