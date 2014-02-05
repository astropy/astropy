from __future__ import print_function
import socket

# save original socket method for restoration
# These are global so that re-calling the turn_off_internet function doesn't
# overwrite them again
socket_original = socket.socket
socket_create_connection = socket.create_connection
socket_bind = socket.socket.bind
socket_connect = socket.socket.connect

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

    socket_functions = {'bind':socket_bind,
                        'connect':socket_connect,
                        'create_connection':socket_create_connection}
    parents = {'bind':socket.socket,
               'connect':socket.socket,
               'create_connection':socket}
    guard_functions = dict([(fnc,blocker_creator(oldfnc, hasattr(parents[fnc],'bind')))
                            for fnc,oldfnc in socket_functions.iteritems()])

    for fnc, guard in guard_functions.iteritems():
        setattr(parents[fnc], fnc, guard)

    return socket

def turn_on_internet(verbose=False):
    """
    Restore internet access.  Not used, but kept in case it is needed.
    """
    if verbose:
        print("Internet access enabled")
    #setattr(socket, 'socket', socket_original)
    setattr(socket, 'create_connection', socket_create_connection)
    setattr(socket.socket, 'bind', socket_bind)
    setattr(socket.socket, 'connect', socket_connect)
    return socket
