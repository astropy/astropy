from __future__ import print_function
import socket

# save original socket method for restoration
socket_original = socket.socket

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
    def guard(*args, **kwargs):
        pytest.fail("An attempt was made to connect to the internet")
    setattr(socket, 'socket', guard)
    return socket

def turn_on_internet(verbose=False):
    """
    Restore internet access.  Not used, but kept in case it is needed.
    """
    if verbose:
        print("Internet access enabled")
    setattr(socket, 'socket', socket_original)
    return socket
