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
            valid_hosts = ('localhost', '127.0.0.1', '::1')
        else:
            host = args[0][0]
            valid_hosts = ('localhost', '127.0.0.1')
        if any([h in host for h in valid_hosts]):
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

    socket.create_connection = check_internet_off(socket_create_connection)
    socket.socket.bind = check_internet_off(socket_bind)
    socket.socket.connect = check_internet_off(socket_connect)

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
