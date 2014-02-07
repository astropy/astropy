from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import socket

# save original socket method for restoration
# These are global so that re-calling the turn_off_internet function doesn't
# overwrite them again
socket_original = socket.socket
socket_create_connection = socket.create_connection
socket_bind = socket.socket.bind
socket_connect = socket.socket.connect


# ::1 is apparently another valid name for localhost?
# it is returned by getaddrinfo when that function is given localhost


def check_internet_off(original_function):
    def new_function(*args, **kwargs):
        if isinstance(args[0], socket.socket):
            host = args[1][0]
        else:
            host = args[0][0]
        if '127.0.0.1' in host or 'localhost' in host or '::1' in host:
            return original_function(*args, **kwargs)
        else:
            raise IOError("An attempt was made to connect to the internet "
                          "by a test that was not marked `remote_data`.")
    return new_function


def turn_off_internet(verbose=False):
    """
    Disable internet access via python by preventing connections from being
    created using the socket module.  Presumably this could be worked around by
    using some other means of accessing the internet, but all default python
    modules (urllib, requests, etc.) use socket [citation needed].
    """
    __tracebackhide__ = True
    if verbose:
        print("Internet access disabled")

    # ::1 is apparently another valid name for localhost?
    # it is returned by getaddrinfo when that function is given localhost

    def blocker_creator(return_fn, isclass=False):
        __tracebackhide__ = True
        def blocker(*args, **kwargs):
            if isclass:
                address = args[1]
            else:
                address = args[0]
            host,port = address[:2]
            if not ('127.0.0.1' in host or
                    'localhost' in host or
                    '::1' in host): 
                raise IOError("An attempt was made to connect to the internet")
            else:
                return return_fn(*args, **kwargs)
        return blocker

    setattr(socket.socket, 'bind', blocker_creator(socket_bind, isclass=True))
    setattr(socket.socket, 'connect', blocker_creator(socket_connect, isclass=True))
    setattr(socket, 'create_connection', blocker_creator(socket_create_connect, isclass=False))

    return socket


def turn_on_internet(verbose=False):
    """
    Restore internet access.  Not used, but kept in case it is needed.
    """
    if verbose:
        print("Internet access enabled")
    socket.create_connection = socket_create_connection
    socket.socket.bind = socket_bind
    socket.socket.connect = socket_connect
    return socket
