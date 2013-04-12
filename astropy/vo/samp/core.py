# Licensed under a 3-clause BSD style license - see LICENSE.rst
#!/usr/bin/env python

###################################################################
##
##   sampy.py - This module contains classes to create a SAMP Hub
##              and interface an application to a running SAMP Hub
##              (Standard Profile)
##
##
##   Copyright (C) 2008  INAF-IASF Milano
##
##   This program is free software; you can redistribute it and/or
##   modify it under the terms of the GNU General Public License
##   as published by the Free Software Foundation; either version 2
##   of the License, or (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
##   MA 02110-1301, USA.
##
##   Authors:
##
##   Luigi Paioro
##   INAF - Istituto Nazionale di Astrofisica
##
##   IASF Milano
##   Via Bassini 15, I-20133 Milano, Italy
##   E-mail: luigi at iasf-milano.inaf.it
##   Site  : http://www.iasf-milano.inaf.it/
##
################################################################################
##
##   Automatic keywords:
##   $Date: 2013-04-08 15:20:13 +0200 (Mon, 08 Apr 2013) $
##   $Revision: 1140 $
##   $Author: luigi $
##   $HeadURL: http://cosmos.lambrate.inaf.it/svn/fase/trunks/sampy/sampy.py $
##   $Id: sampy.py 1140 2013-04-08 13:20:13Z luigi $
##
################################################################################


__author__   = "Luigi Paioro"
__author_email__ = "luigi at iasf-milano.inaf.it"
__url__ = "http://packages.python.org/sampy/"
__description__ = "SAMPy is an IVOA SAMP (Simple Application Messaging Protocol) messaging system implementation in Python."
__long_description__ = """SAMPy is a Python implementation of the `IVOA Simple Application
Messaging Protocol`_ version 1.3.

SAMPy Python module (``sampy.py``) provides classes to easily:

1) instantiate one or multiple Hubs;

2) interface an application or script to a running Hub;

3) create and manage a SAMP client.

SAMPy package provides also a stand-alone program ``sampy`` capable to
instantiate a persistent Hub. In order to have a full description of
``sampy`` stand-alone program options, type the command:

::

  shell > sampy -h


.. _IVOA Simple Application Messaging Protocol: http://www.ivoa.net/Documents/latest/SAMP.html
"""
__status__ = "release"
__release__ = "1.3.0"
__revision__ = "$Revision: 1140 $"
__date__ = "$Date: 2013-04-08 15:20:13 +0200 (Mon, 08 Apr 2013) $"
__profile_version__ = "1.3"

__all__ = ["SAMPHubServer",
           "SAMPHubProxy",
           "SAMPClient",
           "SAMPHubError",
           "SAMPClientError",
           "SAMPProxyError",
           "SAMPLog",
           "SAMPIntegratedClient",
           "SAMPMsgReplierWrapper",
           "SAMP_STATUS_OK",
           "SAMP_STATUS_WARNING",
           "SAMP_STATUS_ERROR",
           "SAMP_HUB_SINGLE_INSTANCE",
           "SAMP_HUB_MULTIPLE_INSTANCE",
           "SAMP_RESTRICT_GROUP",
           "SAMP_RESTRICT_OWNER"]

__doc__ = \
        """
        This module contains classes to create a SAMP Hub and/or interface 
        an application to a running SAMP Hub (Standard Profile).

        SAMPy is very simple to use for the creation of a running SAMP Hub. The only thing
        you have to do is to create a L{SAMPHubServer} instance and start it:

        >>> from sampy import *
        >>> hub = SAMPHubServer()
        >>> hub.start()

        C{sampy.py} module file can also be used as a stand-alone program to launch
        new SAMP Hub instances. Just run it:

        C{shell $: ./sampy.py}

        To have a more exhaustive documentation on C{sampy.py} as stand-alone program,
        please read the C{README} file distributed with SAMPy.

        SAMPy can also be used to easily interface your Python script (or application) to a
        running SAMP Hub and create a callable Client thanks to the L{SAMPIntegratedClient}
        class which integrates the funcionalities of the more specific classes L{SAMPClient}
        and L{SAMPHubProxy}. Here below an example of usage (see also C{iclient.py} and
        C{client.py} files distributed with SAMPy):

        >>> from sampy import *
        >>> import time
        >>>
        >>> # Create a client
        ... cli1 = SAMPIntegratedClient(name = "Client 1", description = "Test Client 1",
        ...                             metadata = {"cli1.version":"0.01"})
        >>> # Create another client (name and description included in metadata argument)
        ... cli2 = SAMPIntegratedClient(metadata = {"samp.name":"Client 2",
        ...                                         "samp.description.text":"Test Client 2",
        ...                                         "cli2.version":"0.25"})
        >>>
        >>> # Connect them
        ... cli1.connect()
        >>> cli2.connect()
        >>>
        >>>
        >>> print "CLI1", cli1.getPrivateKey(), cli1.getPublicId()
        CLI1 0d7f4500225981c104a197c7666a8e4e cli#1
        >>> print "CLI2", cli2.getPrivateKey(), cli2.getPublicId()
        CLI2 72b8ad5ccf0fb3a997d733f6673a960e cli#2
        >>>
        >>>
        >>> # Function called when a notification is received
        ... def test_receive_notification(private_key, sender_id, mtype, params, extra):
        ...   print "Notification:", private_key, sender_id, mtype, params, extra
        ...
        >>> # Function called when a call is received
        ... def test_receive_call(private_key, sender_id, msg_id, mtype, params, extra):
        ...   print "Call:", private_key, sender_id, msg_id, mtype, params, extra
        ...   cli1.ereply(msg_id, SAMP_STATUS_OK, result = {"txt": "printed"})
        ...
        >>> # Function called when a response is received
        ... def test_receive_response(private_key, sender_id, msg_id, response):
        ...   print "Response:", private_key, sender_id, msg_id, response
        ...
        >>> # Subscribe Client 1 to "samp.*" and "samp.app.*" MType and bind it to
        ... # the related functions
        ... cli1.bindReceiveNotification("samp.app.*", test_receive_notification)
        >>> cli1.bindReceiveCall("samp.app.*", test_receive_call)
        >>>
        >>> # Bind Client 2 message-tags received to suitable functions
        ... cli2.bindReceiveResponse("my-dummy-print", test_receive_response)
        >>> cli2.bindReceiveResponse("my-dummy-print-specific", test_receive_response)
        >>> # Client 2 notifies to All "samp.app.echo" MType using myhub
        ... cli2.enotifyAll("samp.app.echo", txt = "Hello world!")
        ['cli#2']
        >>> Notification: 0d7f4500225981c104a197c7666a8e4e cli#2 samp.app.echo {'txt': 'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
        >>>
        >>> print cli2.getSubscribedClients("samp.app.echo")
        {'cli#2': {}}
        >>> # Client 2 calls to All "samp.app.echo" MType using "my-dummy-print"
        ... # as message-tag
        ... print cli2.callAll("my-dummy-print",
        ...                    {"samp.mtype": "samp.app.echo",
        ...                     "samp.params": {"txt": "Hello world!"}})
        {'cli#1': 'msg#1;;cli#hub;;cli#2;;my-dummy-print'}
        >>> Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2 msg#1;;cli#hub;;cli#2;;my-dummy-print samp.app.echo {'txt': 'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
        Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}
        >>>
        >>> # Client 2 calls "samp.app.echo" MType on Client 1 tagging it as
        ... # "my-dummy-print-specific"
        ... try:
        ...   print cli2.call(cli1.getPublicId(),
        ...                   "my-dummy-print-specific",
        ...                   {"samp.mtype": "samp.app.echo",
        ...                    "samp.params": {"txt": "Hello Cli 1!"}})
        ... except SAMPProxyError, e:
        ...   print "Error (%d): %s" % (e.faultCode, e.faultString)
        ...
        msg#2;;cli#hub;;cli#2;;my-dummy-print-specific
        >>> Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2 msg#2;;cli#hub;;cli#2;;my-dummy-print-specific samp.app.echo {'txt': 'Hello Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
        Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print-specific {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}
        >>>
        >>>
        >>>
        >>>
        >>> # Function called to test synchronous calls
        ... def test_receive_sync_call(private_key, sender_id, msg_id, mtype, params, extra):
        ...   print "SYNC Call:", sender_id, msg_id, mtype, params, extra
        ...   time.sleep(2)
        ...   cli1.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
        ...                       "samp.result": {"txt": "printed sync"}})
        ...
        >>>
        >>> # Bind test MType for sync calls
        ... cli1.bindReceiveCall("samp.test", test_receive_sync_call)
        >>>
        >>>
        >>>
        >>>
        >>> try:
        ...   # Sync call
        ...   print cli2.callAndWait(cli1.getPublicId(),
        ...                          {"samp.mtype": "samp.test",
        ...                           "samp.params": {"txt": "Hello SYNCRO Cli 1!"}},
        ...                          "10")
        ... except SAMPProxyError, e:
        ...   # If timeout expires than a SAMPProxyError is returned
        ...   print "Error (%s): %s" % (e.faultCode, e.faultString)
        ...
        SYNC Call: cli#2 msg#3;;cli#hub;;cli#2;;sampy::sync::call samp.test {'txt': 'Hello SYNCRO Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
        {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed sync'}}
        >>>
        >>> cli1.disconnect()
        >>> cli2.disconnect()
        >>>
        >>>
        >>>                               

        """


import os
import sys
import re
import stat
import xmlrpclib
import threading

import datetime, time
import socket
import select
import httplib
import getpass
import hashlib
import base64
import random
import copy
import string
import platform
import inspect
import urllib
import urllib2
import urlparse

import traceback
import StringIO

try:
  import Tkinter
  HAS_TKINTER = True
except:
  HAS_TKINTER = False
  
import Queue


PYTHON_VERSION = float(platform.python_version()[:3])

try:
  import bsddb
except ImportError:
  BDB_SUPPORT = False
else:
  BDB_SUPPORT = True

try:
  import ssl
except ImportError:
  SSL_SUPPORT = False
else:
  SSL_SUPPORT = True

try:
  from urlparse import parse_qs
except ImportError:
  from cgi import parse_qs


from SimpleHTTPServer import *
from SimpleXMLRPCServer import *
from SocketServer import ThreadingMixIn, BaseServer
from BaseHTTPServer import HTTPServer, BaseHTTPRequestHandler

#: General constant for samp.ok status string
SAMP_STATUS_OK = "samp.ok"
#: General constant for samp.warning status string
SAMP_STATUS_WARNING = "samp.warning"
#: General constant for samp.error status string
SAMP_STATUS_ERROR = "samp.error"

#: General constant to specify single instance Hub running mode
SAMP_HUB_SINGLE_INSTANCE = "single"
#: General constant to specify multiple instance Hub running mode
SAMP_HUB_MULTIPLE_INSTANCE = "multiple"

#: General constant to specify the access restriction (through Basic Authentication) to the GROUP
SAMP_RESTRICT_GROUP = "GROUP"
#: General constant to specify the access restriction (through Basic Authentication) to the OWNER
SAMP_RESTRICT_OWNER = "OWNER"


SAMPY_ICON = """iVBORw0KGgoAAAANSUhEUgAAABgAAAAWCAYAAADafVyIAAAABmJLR0QA/wD/AP+gvaeTAAAACXBI
WXMAAA7DAAAOwwHHb6hkAAAAB3RJTUUH2AgFDTkeyUSsmwAAA5FJREFUSMellH9M1GUcx1/PXTfc
MfGGguR1IQeGdZN98ThFSMGmJLIxpFwF0UgXMZXJYm3Qyrm2XLW5wtVGKx1ZOWTgmIgZ5mamaePX
ibsRQQeJ0KEhxK9BePf01133xWOgvv/67nk+3/freT4/HiGlZLZGJkfl0csnaOqxMzo1TvjipZjD
VrIuSiE1Nkks0gWxUInZgN/+6pY5X+7hzthQwB+Cg/Q8t/pZdm98BWtknJi93zXolKuWm0VAgNvj
IaP8VekY6ATAaIhgzRNPo9Pq6B3qo/fvPsamxn3xislCyfOFpMYm+Qzfq/tYFmzKxRRqFACP+dNb
/2z3mW98aj3f7P5MaDUa1QmbeuzyUtc1zjsuYe9zkPdVEYrJIou3FpBoXpvwx51e+kdcmEKN3Afw
mgOkPZPCbHMAW5QibFEKJWmFdA06ZW1LA9XN9eQf2w/QnGi2qlKnAui0uv9r4eqet4CrlptF6fYi
SrcX4RjolACWFbGquqiOmBht9X1XN9Vz3vGTXGi3WFbEitnm9wGiwyJFcowNgBn3DLsq36K09gM5
NDG8YNC8ber657bc9kkOQxPDvrVFuiCy4tPJit9GcoxNPBIAoO9uv3zj67dVRfcqbPFSMpU0MtZs
wRaliIcCAEzNTHO4sUJWXjnJ1Mx0wBijIYL85Jd4eV0WBn2IeCCAf8qOXDhKXds51ZD5y5vCwtTX
iA6LFA8E8GpietJaf72x+fjVGm7c6ggYo9VoyNuwk9L0fQnBQfoWH6Cu7Zysbj7NyOQoANbIOF6w
ZqCYLAGvfeNWhzx+tYa6tu8Dpm/lMhOVu8qJDosUoqb5jCyuOhDwRFstKZSl78P/8fLX0MSwPPxD
BSd+PYXb41HtRSwJ5+z+bxHpn+bKua7svfbOhEwKNuXOCeoadMp3Tn3INWeLaj1vw4uIlI92yJt3
B8hU0vj33gz11xvnhKXGJpGppJEUY8NoiFDB3B4Pe78rkw3tP6pSJcrPfi6P/VzFxXfPYNCHiAsd
l2XJyYOqQZurRVc/HkPEknAA+oddXOz8RRWjmCyIe243rUfel7big8K/NYurDnClu4lH0aHsMjRa
jQZbTAK4XNKvQKLqzQpxKLsMgz7kocwz4raQsz5bzDsHY1PjX5y2NxbUtjbQ1GOf1zg4SM/ezfns
2fy60Go0Cx80L8x+01Hw+6CTrttO2v2678lQI4nmtWTFp6uejf8A94dzfMTZww8AAAAASUVORK5C
YII="""


def internet_on():
  try:
    response = urllib2.urlopen('http://74.125.113.99',timeout=1)
    return True
  except urllib2.URLError as err: pass
  return False

SAFE_MTYPES = ["samp.app.*", "samp.msg.progress", "table.*", "image.*",
               "coord.*", "spectrum.*", "bibcode.*", "voresource.*"]

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
  
  def __init__(self, size, proxy_class, *args, **keywords):
    
    self._proxies = Queue.Queue(size)
    for i in xrange(size):
      self._proxies.put(proxy_class(*args, **keywords))
      
  def __getattr__(self, name):
    # magic method dispatcher
    return _ServerProxyPoolMethod(self._proxies, name)
    

if HAS_TKINTER:
  class WebProfilePopupDialogue(Tkinter.Tk):
  
    def __init__(self, queue, screenName=None, baseName=None, className='Tk',
                 useTk=1, sync=0, use=None):
  
      Tkinter.Tk.__init__(self, screenName, baseName, className,
                          useTk, sync, use)
  
      self._queue = queue
  
      self.title("SAMP Hub")
      self._text = Tkinter.Label(self, font=("Helvetica", 14), \
                                 fg="red", justify=Tkinter.CENTER)
      self._text.pack(padx=5, pady=5, expand=1, fill=Tkinter.X,)
  
      a = Tkinter.Button(self, text="CONSENT", command=self._consent)
      a.pack(padx=5, pady=5, expand=1, fill=Tkinter.X, side=Tkinter.LEFT)
  
      r = Tkinter.Button(self, text="REJECT", command=self._reject)
      r.pack(padx=5, pady=5, expand=1, fill=Tkinter.X, side=Tkinter.RIGHT)
  
      self.protocol("WM_DELETE_WINDOW", self._reject)
  
      self.withdraw()
  
  
    def showPopup(self, request):
      
      samp_name = "unknown"
      
      if isinstance(request[0], str):
        # To support the old protocol version
        samp_name = request[0]
      else:
        samp_name = request[0]["samp.name"]
  
      text = \
  """A Web application which declares to be
  
  Name: %s
  Origin: %s
  
  is requesting to be registered with the SAMP Hub.
  Pay attention that if you permit its registration, such
  application will acquire all current user privileges, like
  file read/write.
  
  Do you give your consent?""" % (samp_name, request[2])
  
      self._text.configure(text=text)
      self.deiconify()
  
    def _consent(self):
      self._queue.put(True)
      self.withdraw()
  
    def _reject(self):
      self._queue.put(False)
      self.withdraw()
  


class SAMPMsgReplierWrapper(object):
  """Decorator class/function that allows to automatically grab
  errors and returned maps (if any) from a function bound
  to a SAMP call (or notify).
  """

  def __init__(self, cli):
    """Decorator initialization, accepting the instance of the
    client that receives the call or notification.

    @param cli: a SAMP client instance.
    @type cli: L{SAMPIntegratedClient} or L{SAMPClient}
    """
    self.cli = cli

  def __call__(self, f):

    def wrapped_f(*args):
      if (inspect.ismethod(f) and f.im_func.func_code.co_argcount == 6) or \
         (inspect.isfunction(f) and f.func_code.co_argcount == 5) or \
         args[2] is None:

        # It is a notification
        f(*args)

      else:
        # It's a call
        try:
          result = f(*args)
          if result:
            self.cli.hub.reply(self.cli.getPrivateKey(), args[2],
                               {"samp.status": SAMP_STATUS_ERROR,
                                "samp.result": result})
        except:
          err = StringIO.StringIO()
          traceback.print_exc(file=err)
          txt = err.getvalue()
          self.cli.hub.reply(self.cli.getPrivateKey(), args[2],
                             {"samp.status": SAMP_STATUS_ERROR,
                              "samp.result": {"txt": txt}})

    return wrapped_f


class SAMPLog(object):
  """
  SAMP Hub logging class. It provides methods for gracefully print SAMPy
  logging messages.
  """

  #: Disable logging at all
  OFF = 0
  #: Log error messages only
  ERROR = 1
  #: Log errors and warnings
  WARNING = 2
  #: Log info, warning and error messages
  INFO = 3 
  #: Log everything for debugging
  DEBUG = 4

  def __init__(self, level = INFO, stdout = sys.stdout, stderr = sys.stderr, prefix = "SAMP"):
    """
    Log class constructor.

    @param level: logging level (one among L{OFF}, L{ERROR}, L{WARNING},
    L{INFO} and L{DEBUG}). By default it is set to L{INFO}.
    @type level: int

    @param stdout: file-like output device. By default it is set to sys.stdout.
    @type stdout: file

    @param stderr: file-like error device. By default it is set to sys.stderr.
    @type stderr: file

    @param prefix: prefix string to logging messages ("SAMP" by default)
    @type prefix: string
    """
    self._level = level
    self._stdout = stdout
    self._stderr = stderr
    self._prefix = prefix


  def setLevel(self, level):
    """
    Set the logging level.

    @param level: one among L{OFF}, L{ERROR}, L{WARNING}, L{INFO} and L{DEBUG}.
    @type level: int
    """
    self._level = level

  def getLevel(self):
    """
    Return the current logging level.

    @return: the current logging level. See L{setLevel}.
    @rtype: int
    """
    return self._level

  def translateLevel(self, level):
    """
    Translate a logging level from the numeric form to a string form and vice versa.
    For example: L{ERROR} is translated in C{"ERROR"} string and vice versa.

    @param level: the logging level to be translated
    @type level: int/str

    @return: the logging level traslated from the numeric form to the string form and vice versa.
    @rtype: int/str
    """

    if isinstance(level, int):
      if level == SAMPLog.OFF:
        return "OFF"
      elif level == SAMPLog.INFO:
        return "INFO"
      elif level == SAMPLog.ERROR:
        return "ERROR"
      elif level == SAMPLog.WARNING:
        return "WARNING"
      elif level == SAMPLog.DEBUG:
        return "DEBUG"
    elif isinstance(level, str):
      if level.upper() == "OFF":
        return SAMPLog.OFF
      elif level.upper() == "INFO":
        return SAMPLog.INFO
      elif level.upper() == "ERROR":
        return SAMPLog.ERROR
      elif level.upper() == "WARNING":
        return SAMPLog.WARNING
      elif level.upper() == "DEBUG":
        return SAMPLog.DEBUG
    else:
      return None


  def info(self, message):
    """
    Log an INFO message.

    @param message: the message to be logged
    @type message: string
    """
    if self._level >= self.INFO:
      self._stdout.write('[%s] Info    (%s): %s\n' % 
                         (self._prefix, datetime.datetime.now().isoformat(), message))
      self._stdout.flush()

  def error(self, message):
    """
    Log an ERROR message.

    @param message: the message to be logged
    @type message: string
    """
    if self._level >= self.ERROR:
      self._stderr.write('[%s] Error   (%s): %s\n' % 
                         (self._prefix, datetime.datetime.now().isoformat(), message))
      self._stderr.flush()

  def warning(self, message):
    """
    Log a WARNING message.

    @param message: the message to be logged
    @type message: string
    """
    if self._level >= self.WARNING:
      self._stderr.write('[%s] Warning (%s): %s\n' % 
                         (self._prefix, datetime.datetime.now().isoformat(), message))
      self._stderr.flush()

  def debug(self, message):
    """
    Log a DEBUG message.

    @param message: the message to be logged
    @type message: string
    """
    if self._level >= self.DEBUG:
      self._stdout.write('[%s] Debug   (%s): %s\n' % 
                         (self._prefix, datetime.datetime.now().isoformat(), message))
      self._stdout.flush()





class SAMPSimpleXMLRPCRequestHandler(SimpleXMLRPCRequestHandler):

  def do_GET(self):

    if self.path == '/sampy/icon':
      self.send_response(200, 'OK')
      self.send_header('Content-Type', 'image/png')
      self.end_headers()
      self.wfile.write(base64.decodestring(SAMPY_ICON))
      
  if PYTHON_VERSION >= 2.7:
    
    def do_POST(self):
        """Handles the HTTP POST request.

        Attempts to interpret all HTTP POST requests as XML-RPC calls,
        which are forwarded to the server's _dispatch method for handling.
        """

        # Check that the path is legal
        if not self.is_rpc_path_valid():
            self.report_404()
            return

        try:
          # Get arguments by reading body of request.
          # We read this in chunks to avoid straining
          # socket.read(); around the 10 or 15Mb mark, some platforms
          # begin to have problems (bug #792570).
          max_chunk_size = 10*1024*1024
          size_remaining = int(self.headers["content-length"])
          L = []
          while size_remaining:
              chunk_size = min(size_remaining, max_chunk_size)
              L.append(self.rfile.read(chunk_size))
              size_remaining -= len(L[-1])
          data = ''.join(L)
          
          params, method = xmlrpclib.loads(data)
          
          if method == "samp.webhub.register":
            params = list(params)
            params.append(self.client_address)
            if self.headers.has_key('Origin'):
              params.append(self.headers.get('Origin'))
            else:
              params.append('unknown')
            params = tuple(params)
            data = xmlrpclib.dumps(params, methodname=method)
  
          elif method in ('samp.hub.notify', 'samp.hub.notifyAll', 
                          'samp.hub.call', 'samp.hub.callAll', 
                          'samp.hub.callAndWait'):
  
            user = "unknown"
  
            if self.headers.has_key('Authorization'):
              # handle Basic authentication
              (enctype, encstr) =  self.headers.get('Authorization').split()
              user, password = base64.standard_b64decode(encstr).split(':')
  
            if method == 'samp.hub.callAndWait':
              params[2]["host"] = self.address_string()
              params[2]["user"] = user
            else:
              params[-1]["host"] = self.address_string()
              params[-1]["user"] = user
  
            data = xmlrpclib.dumps(params, methodname=method)

          data = self.decode_request_content(data)
          if data is None:
              return #response has been sent

          # In previous versions of SimpleXMLRPCServer, _dispatch
          # could be overridden in this class, instead of in
          # SimpleXMLRPCDispatcher. To maintain backwards compatibility,
          # check to see if a subclass implements _dispatch and dispatch
          # using that method if present.
          response = self.server._marshaled_dispatch(
                  data, getattr(self, '_dispatch', None), self.path
              )
        except Exception, e: # This should only happen if the module is buggy
          # internal error, report as HTTP server error
          self.send_response(500)

          # Send information about the exception if requested
          if hasattr(self.server, '_send_traceback_header') and \
                  self.server._send_traceback_header:
              self.send_header("X-exception", str(e))
              self.send_header("X-traceback", traceback.format_exc())

          self.send_header("Content-length", "0")
          self.end_headers()
        else:
          # got a valid XML RPC response
          self.send_response(200)
          self.send_header("Content-type", "text/xml")
          if self.encode_threshold is not None:
              if len(response) > self.encode_threshold:
                  q = self.accept_encodings().get("gzip", 0)
                  if q:
                      try:
                          response = xmlrpclib.gzip_encode(response)
                          self.send_header("Content-Encoding", "gzip")
                      except NotImplementedError:
                          pass
          self.send_header("Content-length", str(len(response)))
          self.end_headers()
          self.wfile.write(response)
    

  elif PYTHON_VERSION >= 2.6 and PYTHON_VERSION < 2.7:

    def do_POST(self):
      """Handles the HTTP POST request.

      Attempts to interpret all HTTP POST requests as XML-RPC calls,
      which are forwarded to the server's _dispatch method for handling.
      """

      # Check that the path is legal
      if not self.is_rpc_path_valid():
        self.report_404()
        return

      try:
        # Get arguments by reading body of request.
        # We read this in chunks to avoid straining
        # socket.read(); around the 10 or 15Mb mark, some platforms
        # begin to have problems (bug #792570).
        max_chunk_size = 10*1024*1024
        size_remaining = int(self.headers["content-length"])
        L = []
        while size_remaining:
          chunk_size = min(size_remaining, max_chunk_size)
          L.append(self.rfile.read(chunk_size))
          size_remaining -= len(L[-1])
        data = ''.join(L)

        params, method = xmlrpclib.loads(data)

        if method == "samp.webhub.register":
          params = list(params)
          params.append(self.client_address)
          if self.headers.has_key('Origin'):
            params.append(self.headers.get('Origin'))
          else:
            params.append('unknown')
          params = tuple(params)
          data = xmlrpclib.dumps(params, methodname=method)

        elif method in ('samp.hub.notify', 'samp.hub.notifyAll', 
                        'samp.hub.call', 'samp.hub.callAll', 
                        'samp.hub.callAndWait'):

          user = "unknown"

          if self.headers.has_key('Authorization'):
            # handle Basic authentication
            (enctype, encstr) =  self.headers.get('Authorization').split()
            user, password = base64.standard_b64decode(encstr).split(':')

          if method == 'samp.hub.callAndWait':
            params[2]["host"] = self.address_string()
            params[2]["user"] = user
          else:
            params[-1]["host"] = self.address_string()
            params[-1]["user"] = user

          data = xmlrpclib.dumps(params, methodname=method)

        # In previous versions of SimpleXMLRPCServer, _dispatch
        # could be overridden in this class, instead of in
        # SimpleXMLRPCDispatcher. To maintain backwards compatibility,
        # check to see if a subclass implements _dispatch and dispatch
        # using that method if present.
        response = self.server._marshaled_dispatch(
          data, getattr(self, '_dispatch', None)
        )
      except Exception, e: # This should only happen if the module is buggy
        # internal error, report as HTTP server error
        self.send_response(500)

        # Send information about the exception if requested
        if hasattr(self.server, '_send_traceback_header') and \
           self.server._send_traceback_header:
          self.send_header("X-exception", str(e))
          self.send_header("X-traceback", traceback.format_exc())

        self.end_headers()
      else:
        # got a valid XML RPC response
        self.send_response(200)
        self.send_header("Content-Type", "text/xml")
        self.send_header("Content-Length", str(len(response)))
        self.end_headers()
        self.wfile.write(response)

        # shut down the connection
        self.wfile.flush()
        self.connection.shutdown(1)

  else:

    def do_POST(self):
      """Handles the HTTP POST request.

      Attempts to interpret all HTTP POST requests as XML-RPC calls,
      which are forwarded to the server's _dispatch method for handling.
      """

      # Check that the path is legal
      if not self.is_rpc_path_valid():
        self.report_404()
        return

      try:
        # Get arguments by reading body of request        self.connection.close().
        # We read this in chunks to avoid straining
        # socket.read(); around the 10 or 15Mb mark, some platforms
        # begin to have problems (bug #792570).
        max_chunk_size = 10*1024*1024
        size_remaining = int(self.headers["content-length"])
        L = []
        while size_remaining:
          chunk_size = min(size_remaining, max_chunk_size)
          L.append(self.rfile.read(chunk_size))
          size_remaining -= len(L[-1])
        data = ''.join(L)

        params, method = xmlrpclib.loads(data)

        if method == "samp.webhub.register":
          params = list(params)
          params.append(self.client_address)
          if self.headers.has_key('Origin'):
            params.append(self.headers.get('Origin'))
          else:
            params.append('unknown')
          params = tuple(params)
          data = xmlrpclib.dumps(params, methodname=method)

        elif method in ('samp.hub.notify', 'samp.hub.notifyAll', 
                        'samp.hub.call', 'samp.hub.callAll', 
                        'samp.hub.callAndWait'):

          user = "unknown"

          if self.headers.has_key('Authorization'):
            # handle Basic authentication
            (enctype, encstr) =  self.headers.get('Authorization').split()
            user, password = base64.standard_b64decode(encstr).split(':')

          if method == 'samp.hub.callAndWait':
            params[2]["host"] = self.address_string()
            params[2]["user"] = user
          else:
            params[-1]["host"] = self.address_string()
            params[-1]["user"] = user

          data = xmlrpclib.dumps(params, methodname=method)

        # In previous versions of SimpleXMLRPCServer, _dispatch
        # could be overridden in this class, instead of in
        # SimpleXMLRPCDispatcher. To maintain backwards compatibility,
        # check to see if a subclass implements _dispatch and dispatch
        # using that method if present.
        response = self.server._marshaled_dispatch(
          data, getattr(self, '_dispatch', None)
        )
      except: # This should only happen if the module is buggy
        # internal error, report as HTTP server error
        self.send_response(500)
        self.end_headers()
      else:
        # got a valid XML RPC response
        self.send_response(200)
        self.send_header("Content-Type", "text/xml")
        self.send_header("Content-Length", str(len(response)))
        self.end_headers()
        self.wfile.write(response)

        # shut down the connection
        self.wfile.flush()
        self.connection.shutdown(1)


class ThreadingXMLRPCServer(ThreadingMixIn, SimpleXMLRPCServer):

  def __init__(self, addr, log = None, requestHandler = SAMPSimpleXMLRPCRequestHandler,
               logRequests = True, allow_none = True, encoding = None):
    self.log = log
    SimpleXMLRPCServer.__init__(self, addr, requestHandler,
                                logRequests, allow_none, encoding)

  def handle_error(self, request, client_address):
    if self.log == None:
      BaseServer.handle_error(self, request, client_address)
    else:
      self.log.warning("Exception happened during processing of request from %s: %s" % (client_address, sys.exc_info()[1]))



class WebProfileRequestHandler(SAMPSimpleXMLRPCRequestHandler):

  def _send_CORS_header(self):

    if not self.headers.get('Origin') is None:

      method = self.headers.get('Access-Control-Request-Method')
      if method and self.command == "OPTIONS":
        # Preflight method
        self.send_header('Content-Length', '0')
        self.send_header('Access-Control-Allow-Origin', self.headers.get('Origin'))
        self.send_header('Access-Control-Allow-Methods', method)
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.send_header('Access-Control-Allow-Credentials', 'true')
      else:
        #Simple method
        self.send_header('Access-Control-Allow-Origin', self.headers.get('Origin'))
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.send_header('Access-Control-Allow-Credentials', 'true')


  def end_headers(self):
    self._send_CORS_header()
    SAMPSimpleXMLRPCRequestHandler.end_headers(self)


  def _serve_cross_domain_xml(self):

    cross_domain = False

    if self.path == "/crossdomain.xml":
      # Adobe standard
      response = """<?xml version='1.0'?>
<!DOCTYPE cross-domain-policy SYSTEM "http://www.adobe.com/xml/dtds/cross-domain-policy.dtd">
<cross-domain-policy>
  <site-control permitted-cross-domain-policies="all"/>
  <allow-access-from domain="*"/>
  <allow-http-request-headers-from domain="*" headers="*"/>
</cross-domain-policy>"""

      self.send_response(200, 'OK')
      self.send_header('Content-Type', 'text/x-cross-domain-policy')
      self.send_header("Content-Length", str(len(response)))
      self.end_headers()
      self.wfile.write(response)
      self.wfile.flush()
      cross_domain = True

    elif self.path == "/clientaccesspolicy.xml":
      # Microsoft standard
      response = """<?xml version='1.0'?>
<access-policy>
  <cross-domain-access>
    <policy>
      <allow-from>
        <domain uri="http://*"/>
      </allow-from>
      <grant-to>
        <resource path="/" include-subpaths="true"/>
      </grant-to>
    </policy>
  </cross-domain-access>
</access-policy>"""

      self.send_response(200, 'OK')
      self.send_header('Content-Type', 'text/xml')
      self.send_header("Content-Length", str(len(response)))
      self.end_headers()
      self.wfile.write(response)
      self.wfile.flush()
      cross_domain = True

    return cross_domain

  def do_POST(self):
    if self._serve_cross_domain_xml():
      return

    return SAMPSimpleXMLRPCRequestHandler.do_POST(self)

  def do_HEAD(self):

    if not self.is_http_path_valid():
      self.report_404()
      return

    if self._serve_cross_domain_xml():
      return


  def do_OPTIONS(self):

    self.send_response(200, 'OK')
    self.end_headers()


  def do_GET(self):

    if not self.is_http_path_valid():
      self.report_404()
      return

    split_path = self.path.split('?')

    if split_path[0] in ['/translator/%s' % clid for clid in self.server.clients]:
      # Request of a file proxying
      urlpath = parse_qs(split_path[1])
      try:
        proxyfile = urllib.urlopen(urlpath["ref"][0])
        self.send_response(200, 'OK')
        self.end_headers()
        self.wfile.write(proxyfile.read())
        proxyfile.close()
      except:
        self.report_404()
        return

    if self._serve_cross_domain_xml():
      return


  def is_http_path_valid(self):

    valid_paths = ["/clientaccesspolicy.xml", "/crossdomain.xml"] + ['/translator/%s' % clid for clid in self.server.clients]
    return self.path.split('?')[0] in valid_paths


class WebProfileXMLRPCServer(ThreadingXMLRPCServer):

  def __init__(self, addr, log = None, requestHandler = WebProfileRequestHandler,
               logRequests = True, allow_none = True, encoding = None):

    self.clients = []
    ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                   logRequests, allow_none, encoding)

  def add_client(self, client_id):
    self.clients.append(client_id)

  def remove_client(self, client_id):
    try:
      self.clients.remove(client_id)
    except:
      pass


if SSL_SUPPORT:

  class HTTPSConnection(httplib.HTTPConnection):
    "This class allows communication via SSL (client side)."

    default_port = httplib.HTTPS_PORT

    def __init__(self, host, port=None, key_file=None, cert_file=None, 
                 cert_reqs=ssl.CERT_NONE, ca_certs=None,
                 ssl_version=ssl.PROTOCOL_SSLv3, strict=None):

      httplib.HTTPConnection.__init__(self, host, port, strict)

      self.key_file = key_file
      self.cert_file = cert_file
      self.cert_reqs = cert_reqs
      self.ca_certs = ca_certs
      self.ssl_version = ssl_version


    def connect(self):
      "Connect to a host on a given (SSL) port."

      sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
      sock.connect((self.host, self.port))
      sslconn = ssl.wrap_socket(sock, server_side = False,
                                certfile = self.cert_file,
                                keyfile = self.key_file,
                                cert_reqs = self.cert_reqs,
                                ca_certs = self.ca_certs,
                                ssl_version = self.ssl_version)
      self.sock = sslconn


  class HTTPS(httplib.HTTP):

    _connection_class = HTTPSConnection

    def __init__(self, host='', port=None, key_file=None, cert_file=None,
                 cert_reqs=ssl.CERT_NONE, ca_certs=None,
                 ssl_version=ssl.PROTOCOL_SSLv3):

      # provide a default host, pass the X509 cert info

      # urf. compensate for bad input.
      if port == 0:
        port = None

      self._setup(self._connection_class(host, port, key_file,
                                         cert_file, cert_reqs,
                                         ca_certs, ssl_version, None))

      # we never actually use these for anything, but we keep them
      # here for compatibility with post-1.5.2 CVS.
      self.key_file = key_file
      self.cert_file = cert_file
      
    def getresponse(self, buffering=False):
      "Get the response from the server."
      return self._conn.getresponse(buffering)


  class SafeTransport(xmlrpclib.Transport):
    """Handles an HTTPS transaction to an XML-RPC server."""

    def __init__(self, key_file=None, cert_file=None,
                 cert_reqs=ssl.CERT_NONE, ca_certs=None,
                 ssl_version=ssl.PROTOCOL_SSLv3, strict=None,
                 use_datetime=0):

      xmlrpclib.Transport.__init__(self, use_datetime)
      self.key_file = key_file
      self.cert_file = cert_file
      self.cert_reqs = cert_reqs
      self.ca_certs = ca_certs
      self.ssl_version = ssl_version


    def make_connection(self, host):
      # create a HTTPS connection object from a host descriptor
      # host may be a string, or a (host, x509-dict) tuple
      host, extra_headers, x509 = self.get_host_info(host)
      return HTTPS(host, None, self.key_file, self.cert_file,
                   self.cert_reqs, self.ca_certs, self.ssl_version)


  class SecureXMLRPCServer(ThreadingXMLRPCServer):

    def __init__(self, addr, keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                 log = None, requestHandler = SimpleXMLRPCRequestHandler,
                 logRequests = True, allow_none = True, encoding = None):
      """
      Secure XML-RPC server.

      It it very similar to SimpleXMLRPCServer but it uses HTTPS for transporting XML data.
      """
      self.keyfile  = keyfile
      self.certfile = certfile
      self.cert_reqs = cert_reqs
      self.ca_certs = ca_certs
      self.ssl_version = ssl_version
      self.allow_reuse_address = True

      ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                     logRequests, allow_none, encoding)

    def get_request (self):
      # override this to wrap socket with SSL
      sock, addr = self.socket.accept()
      sslconn = ssl.wrap_socket(sock, server_side = True,
                                certfile = self.certfile,
                                keyfile = self.keyfile,
                                cert_reqs = self.cert_reqs,
                                ca_certs = self.ca_certs,
                                ssl_version = self.ssl_version)
      return sslconn, addr



if BDB_SUPPORT:

  class BasicAuthSimpleXMLRPCRequestHandler(SAMPSimpleXMLRPCRequestHandler):
    """XML-RPC Request Handler for Basic Authentication support."""

    def __init__(self, request, client_address, server, auth_file, access_restrict = None):
      """
      Constructor.

      @param auth_file: Authentication file path. It is a Berkeley DB file in Hash
      format containing	a set of key=value pairs of the form:
      C{<user name>=md5(<password>)<group 1>,<group 2>,<group 3>,...}.
      @type auth_file: string

      @param access_restrict: Dictionary containing the restriction rules for authentication.
      If the access must be restricted to a specific user then C{access_restrict} is a dictionary
      containing C{{"user"; <user name>}}. If the access must be restricted to the 
      users belonging to a certain group, the C{access_restrict} is a dictionary containing
      C{{"group"; <group name>}}. An additional key can be present: C{"admin": <administrator user>}.
      It defines the name of the administrator user with full access permission.
      @type access_restrict: dictionary
      """

      self.db = bsddb.hashopen(auth_file, "r")
      self.access_restrict = access_restrict
      SimpleXMLRPCRequestHandler.__init__(self, request, client_address, server)
      self.db.close()

    def checkId(self, id, pwd):

      if id in self.db.keys():

        pwdhash = self.db[id][0:16]
        groups =  self.db[id][16:]
        pwd = hashlib.md5(pwd).digest()

        if self.access_restrict != None:

          # ADMIN TEST
          if self.access_restrict.has_key("admin"):
            admin = self.access_restrict["admin"]
            if self.db.has_key(admin):
              adminpwdhash = self.db[admin][0:16]
              if admin == id and adminpwdhash == pwd:
                return True

          # TEST USER RESTRICTION
          if self.access_restrict.has_key("user"):
            if self.access_restrict["user"] == id and pwdhash == pwd:
              return True
            else:
              return False

          # TEST GROUP RESTRICTION
          if self.access_restrict.has_key("group"):
            if self.access_restrict["group"] in groups.split(",") and pwdhash == pwd:
              return True
            else:
              return False
        else:
          if pwdhash == pwd:
            return True
          else:
            return False
      else:
        return False

    def authenticate_client(self):
      validuser = False

      if self.headers.has_key('Authorization'):
        # handle Basic authentication
        (enctype, encstr) =  self.headers.get('Authorization').split()
        (user, password) = base64.standard_b64decode(encstr).split(':')
        validuser = self.checkId(user, password)

      return validuser

    def do_POST(self):    

      if self.authenticate_client():
        SAMPSimpleXMLRPCRequestHandler.do_POST(self)
      else:
        self.report_401()

    def report_401(self):
      # Report a 401 error
      self.send_response(401)
      self.send_header("WWW-Authenticate", "Basic realm=\"Protected access\"");
      self.end_headers()
      # shut down the connection
      self.connection.shutdown(1)
      self.connection.close()


  class BasicAuthXMLRPCServer(ThreadingXMLRPCServer):
    """XML-RPC server with Basic Authentication support."""

    def __init__(self, addr, auth_file, access_restrict = None, log = None,
                 requestHandler = BasicAuthSimpleXMLRPCRequestHandler,
                 logRequests = True, allow_none = True, encoding = None):
      """
      Constructor.

      @param auth_file: Authentication file path. It is a Berkeley DB file in Hash
      format containing	a set of key=value pairs of the form:
      C{<user name>=md5(<password>)<group 1>,<group 2>,<group 3>,...}.
      @type auth_file: string

      @param access_restrict: Dictionary containing the restriction rules for authentication.
      If the access must be restricted to a specific user then access_restrict is a dictionary
      containing C{{"user"; <user name>}}. If the access must be restricted to the 
      users belonging to a certain group, the access_restrict is a dictionary containing
      C{{"group"; <group name>}}. An additional key can be present: C{"admin": <administrator user>}.
      It defines the name of the administrator user with full access permission.
      @type access_restrict: dictionary
      """

      self.auth_file = auth_file
      self.access_restrict = access_restrict

      ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                     logRequests, allow_none, encoding)

    def finish_request(self, request, client_address):
      if self.auth_file != None and self.RequestHandlerClass == BasicAuthSimpleXMLRPCRequestHandler:
        self.RequestHandlerClass(request, client_address, self,
                                 self.auth_file, self.access_restrict)
      else:
        ThreadingXMLRPCServer.finish_request(self, request, client_address)

if SSL_SUPPORT and BDB_SUPPORT:

  class BasicAuthSecureXMLRPCServer(ThreadingXMLRPCServer):
    def __init__(self, addr, keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                 auth_file, access_restrict = None, log = None,
                 requestHandler = BasicAuthSimpleXMLRPCRequestHandler,
                 logRequests = True, allow_none = True, encoding = None):

      self.keyfile  = keyfile
      self.certfile = certfile
      self.cert_reqs = cert_reqs
      self.ca_certs = ca_certs
      self.ssl_version = ssl_version
      self.allow_reuse_address = True
      self.auth_file = auth_file
      self.access_restrict = access_restrict

      ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                     logRequests, allow_none, encoding)

    def get_request (self):
      # override this to wrap socket with SSL
      sock, addr = self.socket.accept()
      sslconn = ssl.wrap_socket(sock, server_side = True,
                                certfile = self.certfile,
                                keyfile = self.keyfile,
                                cert_reqs = self.cert_reqs,
                                ca_certs = self.ca_certs,
                                ssl_version = self.ssl_version)
      return sslconn, addr

    def finish_request(self, request, client_address):
      if self.auth_file != None and self.RequestHandlerClass == BasicAuthSimpleXMLRPCRequestHandler:
        self.RequestHandlerClass(request, client_address, self,
                                 self.auth_file, self.access_restrict)
      else:
        ThreadingXMLRPCServer.finish_request(self, request, client_address)



class SAMPHubServer(object):
  """
  SAMP Hub Server implementation (Standard Profile v1.0)
  """

  def __init__(self, secret = None, addr=None, port=0, lockfile=None, timeout = 0, \
               client_timeout = 0, log = None,
               mode = SAMP_HUB_SINGLE_INSTANCE, label = "",
               owner = "", owner_group = "", auth_file = None,
               access_restrict = None, admin = "admin", https = False,
               keyfile = None, certfile = None,
               cert_reqs = 0, ca_certs = None, ssl_version = 2, web_profile = True,
               pool_size=20):
    """
    SAMP Hub server constructor. The SSL parameters are usable only if Python has
    been compiled with SSL support and/or U{ssl <http://docs.python.org/dev/library/ssl.html>}
    module is installed (available by default since Python 2.6).

    @param secret: Secret code.
    @type secret: string

    @param addr: Listening address (or IP)
    @type addr: string

    @param port: Listening port number.
    @type port: int

    @param lockfile: Custom lockfile name.
    @type port: string

    @param timeout: Hub inactivity timeout. If C{timeout} > 0 then the Hub automatically
    stops after an inactivity period longer than C{timeout} seconds. By default C{timeout}
    is set to 0 (Hub never expires).
    @type timeout: int

    @param client_timeout: Client inactivity timeout. If C{client_timeout} > 0 then the
    Hub automatically	unregisters the clients which result inactive for a period longer
    than C{client_timeout} seconds. By default C{client_timeout} is set to 0 (clients never
    expire).
    @type client_timeout: int

    @param log: a L{SAMPLog} instance for the Hub logging. If missing then
    a standard L{SAMPLog} instance is used.
    @type log: L{SAMPLog}

    @param mode: Defines the Hub running mode. If C{mode} is 'single' then the Hub runs
    using the standard {.samp} lock-file, having a single instance for user desktop
    session. Otherwise, if C{mode} is 'multiple', then the Hub runs using a non-standard
    lock-file, placed in C{.samp-1} directory, of the form C{samp-hub-<PID>-<ID>}, where
    C{<PID>} is the process id and C{<ID>} is a general sub-id (integer number).
    @type mode: string

    @param label: A string used to label the Hub with a human readable name. This string
    is written in the lock-file assigned to the C{hub.label} token.
    @type label: string

    @param owner: General purpose Hub owner name. This value is written in the lock-file
    and assigned to the C{hub.owner.name} token.
    @type owner: string

    @param owner_group: General purpose Hub owner group name. This value is written in the
    lock-file	and assigned to the C{hub.owner.group} token.
    @type owner_group: string

    @param auth_file: Authentication file path used for Basic Authentication. The authentication file
    must be a Berkeley DB file in Hash format containing a set of C{<user name>=md5(<password>)<group
    1>,<group 2>,<group 3>,...)} key/value pairs.
    @type auth_file: string

    @param access_restrict: Define whether the Hub access must be restricted to the Hub owner, to a certain owner
    group or not restricted at all. Values accepted: L{SAMP_RESTRICT_OWNER}, L{SAMP_RESTRICT_GROUP}, B{None}
    @type access_restrict: string

    @param admin: Define the name of the administrator user in case of restricted access. The administrator user
    can always access the hub instance even if it is running with L{SAMP_RESTRICT_OWNER} policy. The administrator
    must be declared in the authentication file.
    @type admin: string

    @param https: set the Hub running on a Secure Sockets Layer connection (HTTPS).
    By default SSL is desabled.
    @type https: boolean

    @param keyfile: Set the file containing the private key for SSL connections. If the
    certificate file (C{certfile}) contains the private key, then C{keyfile} can be omitted.
    @type keyfile: string

    @param certfile: Specify the file which contains a certificate to be used to identify the
    local side of the secure connection.
    @type certfile: string

    @param cert_reqs: The parameter C{cert_reqs} specifies whether a certificate is required
    from the client side of the connection, and whether it will be validated if provided. It
    must be one of the three values L{ssl.CERT_NONE} (certificates ignored), L{ssl.CERT_OPTIONAL}
    (not required, but validated if provided), or L{ssl.CERT_REQUIRED} (required and validated).
    If the value of this parameter is not L{ssl.CERT_NONE}, then the C{ca_certs} parameter must
    point to a file of CA certificates.
    @type cert_reqs: int

    @param ca_certs: The C{ca_certs} file contains a set of concatenated \"Certification Authority\" 
    certificates, which are used to validate certificates passed from the client end of the 
    connection.
    @type ca_certs: string

    @param ssl_version: The C{ssl_version} option specifies which version of the SSL protocol to use.
    Typically, the server chooses	a particular protocol version, and the client must adapt to the
    server's choice. Most of the versions are	not interoperable with the other versions. If not
    specified the default SSL version is  L{ssl.PROTOCOL_SSLv23}. This version provides the most
    compatibility with other versions client side. Other SSL protocol versions are:                        
    L{ssl.PROTOCOL_SSLv2}, L{ssl.PROTOCOL_SSLv3} and L{ssl.PROTOCOL_TLSv1}.
    @type ssl_version: int

    @param web_profile: The C{web_profile} option enables/disables the Web Profile support.
    @type ssl_version: boolean
    
    @param pool_size: The number of socket connections opened to communicate with the clients.
    @type pool_size: int
    """

    # General settings
    self._lockfilename = lockfile
    self._admin = admin
    self._addr = addr
    self._port = port
    self._mode = mode
    self._label = label
    self._owner = owner
    self._owner_group = owner_group
    self._timeout = timeout
    self._client_timeout = client_timeout
    if log == None:
      self._log = SAMPLog()
    else:
      self._log = log
    self._pool_size = pool_size

    # Web Profile variables
    if not HAS_TKINTER:
      if web_profile:
        self._log.info("Tkinter python module is missing. Impossible to start the Web profile.")
        web_profile = False
    
    
    self._web_profile = web_profile
    self._web_profile_server = None
    self._web_profile_callbacks = {}
    self._web_profile_popup_dialogue = None
    self._web_profile_requests_queue = Queue.Queue(1)
    self._web_profile_requests_result = Queue.Queue(1)
    self._web_profile_requests_semaphore = Queue.Queue(1)
    if web_profile:
      try:
        self._web_profile_server = WebProfileXMLRPCServer(('localhost', 21012), self._log, \
                                                          logRequests = False, allow_none = True)
        self._web_profile_server.register_introspection_functions()
        self._log.info("Hub set to run with Web Profile support enabled.")
      except:
        self._log.error("Port 21012 already in use. Impossible to run the Hub with Web Profile support.")
        self._web_profile = web_profile = False

    # SSL general settings
    self._https = https
    self._keyfile = keyfile
    self._certfile = certfile
    self._cert_reqs = cert_reqs
    self._ca_certs = cert_reqs
    self._ssl_version = ssl_version
    # Basic Authentication settings
    self._auth_file = auth_file
    self._access_restrict = access_restrict

    # Reformat access_restrict string to suitable dictionary
    if access_restrict != None:
      if access_restrict == SAMP_RESTRICT_GROUP:
        access_restrict = {"group": owner_group, "admin": admin}
      elif access_restrict == SAMP_RESTRICT_OWNER:
        access_restrict = {"user": owner, "admin": admin}
      else:
        access_restrict = None

    # Athentication file test
    if auth_file != None:
      if not os.path.isfile(auth_file):
        self._log.error("Unable to load authentication file!")
        raise SAMPHubError("Unable to load authentication file!")


    self._host_name = "127.0.0.1"
    if internet_on():
      try:
        self._host_name = socket.getfqdn()
        socket.getaddrinfo(self._addr or self._host_name, self._port or 0)
      except:
        pass

    # XML-RPC server settings
    if https:

      if keyfile != None and not os.path.isfile(keyfile):
        self._log.error("Unable to load SSL private key file!")
        raise SAMPHubError("Unable to load SSL private key file!")

      if certfile == None or not os.path.isfile(certfile):
        self._log.error("Unable to load SSL cert file!")
        raise SAMPHubError("Unable to load SSL cert file!")

      if auth_file != None:
        self._log.info("Hub set for Basic Authentication using SSL.")
        self._server = BasicAuthSecureXMLRPCServer((self._addr or self._host_name, self._port or 0), 
                                                   keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                                                   auth_file, access_restrict, self._log,
                                                   logRequests = False, allow_none = True)
      else:
        self._log.info("Hub set for using SSL.")
        self._server = SecureXMLRPCServer((self._addr or self._host_name, self._port or 0), 
                                          keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                                          self._log, logRequests = False, allow_none = True)

      self._port = self._server.socket.getsockname()[1]
      self._url = "https://%s:%s" %(self._addr or self._host_name,
                                    self._port)
    else:

      if auth_file != None:
        self._log.info("Hub set for Basic Authentication.")
        self._server = BasicAuthXMLRPCServer((self._addr or self._host_name, self._port or 0),
                                             auth_file, access_restrict, self._log,
                                             logRequests = False, allow_none = True)
      else:
        self._server = ThreadingXMLRPCServer((self._addr or self._host_name, self._port or 0),
                                             self._log, logRequests = False, allow_none = True)


      self._port = self._server.socket.getsockname()[1]
      self._url = "http://%s:%s" %(self._addr or self._host_name,
                                   self._port)

    self._server.register_introspection_functions()


    # Threading stuff
    self._thread_lock = threading.Lock()
    self._thread_run = None
    self._thread_hub_timeout = None
    self._thread_client_timeout = None
    self._is_running = False

    # Variables for timeout testing:
    self._last_activity_time = None
    self._client_activity_time = {}


    # Hub message id counter, used to create hub msg ids
    self._hub_msg_id_counter = 0  
    # Hub secred code
    self._hub_secret_code_customized = secret
    self._hub_secret = self._createSecretCode()
    # Hub public id (as SAMP client)
    self._hub_public_id = ""
    # Client ids
    # {private_key: (public_id, timestamp)}
    self._private_keys = {}
    # Metadata per client
    # {private_key: metadata}
    self._metadata = {}
    # List of subscribed clients per MType
    # {mtype: private_key list}
    self._mtype2ids = {}
    # List of subscribed MTypes per client
    # {private_key: mtype list}
    self._id2mtypes = {}
    # List of XML-RPC addresses per client
    # {public_id: (XML-RPC address, ServerProxyPool instance)}
    self._xmlrpcEndpoints = {}
    # Synchronous message id heap
    self._sync_msg_ids_heap = {}
    # Public ids counter
    self._client_id_counter = -1

    # Standard Profile only operations
    self._server.register_function(self._ping, 'samp.hub.ping')
    self._server.register_function(self._setXmlrpcCallback, 'samp.hub.setXmlrpcCallback')
    # Stadard API operations
    self._server.register_function(self._register, 'samp.hub.register')
    self._server.register_function(self._unregister, 'samp.hub.unregister')
    self._server.register_function(self._declareMetadata, 'samp.hub.declareMetadata')
    self._server.register_function(self._getMetadata, 'samp.hub.getMetadata')
    self._server.register_function(self._declareSubscriptions, 'samp.hub.declareSubscriptions')
    self._server.register_function(self._getSubscriptions, 'samp.hub.getSubscriptions')
    self._server.register_function(self._getRegisteredClients, 'samp.hub.getRegisteredClients')
    self._server.register_function(self._getSubscribedClients, 'samp.hub.getSubscribedClients')
    self._server.register_function(self._notify, 'samp.hub.notify')
    self._server.register_function(self._notifyAll, 'samp.hub.notifyAll')
    self._server.register_function(self._call, 'samp.hub.call')
    self._server.register_function(self._callAll, 'samp.hub.callAll')
    self._server.register_function(self._callAndWait, 'samp.hub.callAndWait')
    self._server.register_function(self._reply, 'samp.hub.reply')
    # Hub as client operations (see _hubAsClientRequestHandler)
    #self._server.register_function(self._receiveNotification, 'samp.client.receiveNotification')
    #self._server.register_function(self._receiveCall, 'samp.client.receiveCall')
    #self._server.register_function(self._receiveResponse, 'samp.client.receiveResponse')

    if web_profile:
      # Web Profile methods like Standard Profile
      self._web_profile_server.register_function(self._ping, 'samp.webhub.ping')
      self._web_profile_server.register_function(self._unregister, 'samp.webhub.unregister')
      self._web_profile_server.register_function(self._declareMetadata, 'samp.webhub.declareMetadata')
      self._web_profile_server.register_function(self._getMetadata, 'samp.webhub.getMetadata')
      self._web_profile_server.register_function(self._declareSubscriptions, 'samp.webhub.declareSubscriptions')
      self._web_profile_server.register_function(self._getSubscriptions, 'samp.webhub.getSubscriptions')
      self._web_profile_server.register_function(self._getRegisteredClients, 'samp.webhub.getRegisteredClients')
      self._web_profile_server.register_function(self._getSubscribedClients, 'samp.webhub.getSubscribedClients')
      self._web_profile_server.register_function(self._notify, 'samp.webhub.notify')
      self._web_profile_server.register_function(self._notifyAll, 'samp.webhub.notifyAll')
      self._web_profile_server.register_function(self._call, 'samp.webhub.call')
      self._web_profile_server.register_function(self._callAll, 'samp.webhub.callAll')
      self._web_profile_server.register_function(self._callAndWait, 'samp.webhub.callAndWait')
      self._web_profile_server.register_function(self._reply, 'samp.webhub.reply')
      # Methods peculiar for Web Profile
      self._web_profile_server.register_function(self._web_profile_register, 'samp.webhub.register')      
      self._web_profile_server.register_function(self._web_profile_allowReverseCallbacks, 'samp.webhub.allowReverseCallbacks')
      self._web_profile_server.register_function(self._web_profile_pullCallbacks, 'samp.webhub.pullCallbacks')


  def __del__(self):
    self.stop()

  def _timeoutTestHub(self):

    while self._is_running:
      time.sleep(1)
      self._thread_lock.acquire()
      if self._timeout > 0 and self._last_activity_time != None:
        if time.time() - self._last_activity_time >= self._timeout:
          self._thread_lock.release()
          self._log.warning("Timeout expired, Hub is shutting down!")
          self.stop()
          break
      if self._thread_lock.locked() == True:
        self._thread_lock.release()

  def _timeoutTestClient(self):

    while self._is_running:
      time.sleep(1)
      if self._client_timeout > 0:
        now = time.time()
        for private_key in self._client_activity_time.keys():
          if now - self._client_activity_time[private_key] > self._client_timeout \
             and private_key != self._hub_private_key:
            self._log.warning("Client %s timeout expired!" % private_key)
            self._notifyDisconnection(private_key)
            self._unregister(private_key)


  def _hubAsClientRequestHandler(self, method, args):
    if method == 'samp.client.receiveCall':
      return self._receiveCall(*args)
    elif method == 'samp.client.receiveNotification':
      return self._receiveNotification(*args)
    elif method == 'samp.client.receiveResponse':
      return self._receiveResponse(*args)
    elif method == 'samp.app.ping':
      return self._ping(*args)
    else:
      return _hubAsClientRequestHandler

  def _setupHubAsClient(self):
    result = self._register(self._hub_secret)
    self._hub_public_id = result["samp.self-id"]
    self._hub_private_key = result["samp.private-key"]
    self._setXmlrpcCallback(self._hub_private_key, self._url)
    self._declareMetadata(self._hub_private_key, {"samp.name": "Hub",
                                                  "samp.description.text": self._label,
                                                  "author.name": "Luigi Paioro",
                                                  "author.email": "luigi@iasf-milano.inaf.it",
                                                  "author.affiliation": "INAF-IASF Milano",
                                                  "samp.documentation.url": "http://packages.python.org/sampy/",
                                                  "samp.icon.url": self._url + "/sampy/icon"})
    self._declareSubscriptions(self._hub_private_key, {"samp.app.ping":{}, "x-samp.query.by-meta":{}})

  def start(self, wait = False):
    """
    Start the current SAMP Hub instance and create the lock file. Hub start-up can
    be blocking or non blocking depending on the C{wait} parameter.

    @param wait: if True then the Hub process is joined with the caller, blocking the
    code flow. Usually True option is used to run a stand-alone Hub in
    an executable script. If False (default), then the Hub process runs in a
    separated thread. False is usually used in a Python shell.
    @type wait: boolean
    """

    if self._is_running == False:

      self._is_running = True
      self._updateLastActivityTime()

      if self._createLockFile() == False:
        self._is_running = False
        return

      self._setupHubAsClient()
      self._startThreads()

      self._log.info("Hub started")

    if wait and self._is_running:
      self._thread_run.join()

  def _createLockFile(self):

    # Remove lock-files of dead hubs
    self._removeGarbageLockFiles()

    lockfilename = ""
    lockfiledir = ""


    # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
    if os.environ.has_key("SAMP_HUB"):
      # For the time being I assume just the std profile supported.
      if os.environ["SAMP_HUB"].startswith("std-lockurl:"):

        lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:"):]
        lockfile_parsed = urlparse.urlparse(lockfilename)

        if lockfile_parsed[0] != 'file':
          self._log.warning("Unable to start a Hub with lockfile %s. Start-up process aborted." % lockfilename)
          return False
        else:
          lockfilename = lockfile_parsed[2]
    else:
      # If it is a fresh Hub instance
      if self._lockfilename is None:

        self._log.debug("Running mode: " + self._mode)

        if self._mode == SAMP_HUB_SINGLE_INSTANCE:
          lockfilename = ".samp"
        else:
          lockfilename = "samp-hub-%d-%s" % (os.getpid(), threading._counter + 1)

        if os.environ.has_key("HOME"):
          # UNIX
          lockfiledir = os.environ["HOME"]
        else:
          # Windows
          lockfiledir = os.environ["USERPROFILE"]

        if self._mode == SAMP_HUB_MULTIPLE_INSTANCE:
          lockfiledir = os.path.join(lockfiledir, ".samp-1")

        # If missing create .samp-1 directory
        if not os.path.isdir(lockfiledir):
          os.mkdir(lockfiledir)
          os.chmod(lockfiledir, stat.S_IREAD + stat.S_IWRITE + stat.S_IEXEC)

        lockfilename = os.path.join(lockfiledir, lockfilename)

      else:
        self._log.debug("Running mode: multiple")
        lockfilename = self._lockfilename


    hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)

    if hub_is_running:
      self._log.warning("Another SAMP Hub is already running. Start-up process aborted.")
      return False

    self._log.debug("Lock-file: " + lockfilename)

    result = self._new_lockfile(lockfilename)
    if result:
      self._lockfilename = lockfilename

    return result


  def _new_lockfile(self, lockfilename):

    lockfile = open(lockfilename, "w")
    lockfile.close()
    os.chmod(lockfilename, stat.S_IREAD + stat.S_IWRITE)
    lockfile = open(lockfilename, "w")
    lockfile.write("# SAMP lockfile written on %s\n" % datetime.datetime.now().isoformat())
    lockfile.write("# Standard Profile required keys\n")
    lockfile.write("samp.secret=%s\n" % self._hub_secret)
    lockfile.write("samp.hub.xmlrpc.url=%s\n" % self._url)
    lockfile.write("samp.profile.version=%s\n" % __profile_version__)

    # Custom tokens

    lockfile.write("hub.id=%d-%s\n" % (os.getpid(), _THREAD_STARTED_COUNT))

    if self._label == "":
      self._label = "Hub %d-%s" % (os.getpid(), _THREAD_STARTED_COUNT)
    if self._label != "":
      lockfile.write("hub.label=%s\n" % self._label)
    if self._owner != "":	
      lockfile.write("hub.owner.name=%s\n" % self._owner)
    if self._owner_group != "":	
      lockfile.write("hub.owner.group=%s\n" % self._owner_group)

    if self._auth_file != None:
      lockfile.write("hub.access.auth.file=%s\n" % self._auth_file)

    if self._access_restrict != None:
      lockfile.write("hub.access.auth.restrict=%s\n" % self._access_restrict)

    if SSL_SUPPORT and self._https:
      # Certificate request
      cert_reqs_types = ["NONE", "OPTIONAL", "REQUIRED"]
      lockfile.write("hub.ssl.certificate=%s\n" % cert_reqs_types[self._cert_reqs])
      # SSL protocol version
      ssl_protocol_types = ["SSLv2", "SSLv3", "SSLv23", "TLSv1"]
      lockfile.write("hub.ssl.protocol=%s\n" % ssl_protocol_types[self._ssl_version])

    lockfile.close()

    return True


  def _removeGarbageLockFiles(self):

    lockfilename = ""

    # HUB SINGLE INSTANCE MODE

    if os.environ.has_key("HOME"):
      # UNIX
      lockfilename = os.path.join(os.environ["HOME"], ".samp")
    else:
      # Windows
      lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")

    hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)

    if not hub_is_running:
      # If lockfilename belongs to a dead hub, then it is deleted
      if os.path.isfile(lockfilename):
        try:
          os.remove(lockfilename)
        except:
          pass


    # HUB MULTIPLE INSTANCE MODE

    lockfiledir = ""

    if os.environ.has_key("HOME"):
      # UNIX
      lockfiledir = os.path.join(os.environ["HOME"], ".samp-1")
    else:
      # Windows
      lockfiledir = os.path.join(os.environ["USERPROFILE"], ".samp-1")

    if os.path.isdir(lockfiledir):
      for filename in os.listdir(lockfiledir):
        if re.match('samp\\-hub\\-\d+\\-\d+', filename) != None:
          lockfilename = os.path.join(lockfiledir, filename)
          hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)
          if not hub_is_running:
            # If lockfilename belongs to a dead hub, then it is deleted
            if os.path.isfile(lockfilename):
              try:
                os.remove(lockfilename)
              except:
                pass

  def _startThreads(self):
    self._thread_run = threading.Thread(target = self._serve_forever)
    self._thread_run.setDaemon(True)
    self._thread_hub_timeout = threading.Thread(target = self._timeoutTestHub,
                                                name = "Hub timeout test")
    self._thread_client_timeout = threading.Thread(target = self._timeoutTestClient,
                                                   name = "Client timeout test")
    self._thread_run.start()
    self._thread_hub_timeout.start()
    self._thread_client_timeout.start()


  @staticmethod
  def checkRunningHub(lockfilename):
    """
    Test whether a Hub identified by C{lockfilename} is running or not.

    @param lockfilename: lock-file name (path + file name) of the Hub to be tested.
    @type lockfilename: string
    """

    is_running = False
    lockfiledict = {}

    # Check whether a lockfile alredy exists
    try:
      lockfile = urllib.urlopen(lockfilename)
      lockfile_content = lockfile.readlines()
      lockfile.close()
    except:
      return is_running, lockfiledict

    for line in lockfile_content:
      if line.strip()[0] != "#":
        kw, val = line.split("=")
        lockfiledict[kw.strip()] = val.strip()

    if lockfiledict.has_key("samp.hub.xmlrpc.url"):
      try:
        proxy = xmlrpclib.ServerProxy(lockfiledict["samp.hub.xmlrpc.url"].replace("\\", ""),
                                      allow_none=1)
        proxy.samp.hub.ping()
        is_running = True
      except xmlrpclib.ProtocolError:
        # There is a protocol error (e.g. for authentication required),
        # but the server is alive
        is_running = True
      except:
        if SSL_SUPPORT:
          if sys.exc_info()[0] == ssl.SSLError:
            # SSL connection refused for certifcate reasons...
            # anyway the server is alive
            is_running = True						

    return is_running, lockfiledict


  def _createSecretCode(self):
    if self._hub_secret_code_customized != None:
      return self._hub_secret_code_customized
    else:
      now = datetime.datetime.utcnow().isoformat()
      host = self._host_name
      user = getpass.getuser()
      return hashlib.sha1(host + user + now).hexdigest()


  def stop(self):
    """
    Stop the current SAMP Hub instance and delete the lock file.
    """

    if self._is_running:

      self._log.info("Hub is stopping...")

      self._notifyShutdown()

      self._is_running = False

      if (os.path.isfile(self._lockfilename)):
        lockfile = open(self._lockfilename, "r")
        lockfile_content = lockfile.readlines()
        lockfile.close()
        for line in lockfile_content:
          if line.strip()[0] != "#":
            kw, val = line.split("=")
            if kw.strip() == "samp.secret" and val.strip() == self._hub_secret:
              os.remove(self._lockfilename)
              break

    # Reset vaiables
    self._joinAllThreads()

    self._hub_msg_id_counter = 0
    self._hub_secret = self._createSecretCode()
    self._hub_public_id = ""
    self._metadata = {}
    self._private_keys = {}
    self._mtype2ids = {}
    self._id2mtypes = {}
    self._xmlrpcEndpoints = {}
    self._last_activity_time = None
    self._lockfilename = None

    self._log.info("Hub stopped.")

  def _joinOneThread(self, thread_name, timeout):
    t = getattr(self, thread_name)
    if t is None:
      return
    t.join(timeout)
    setattr(self, thread_name, None)

  def _joinAllThreads(self, timeout=1):
    for thread_name in [
      "_thread_run", 
      "_thread_hub_timeout",
      "_thread_client_timeout"]:
      self._joinOneThread(thread_name, timeout)

  def isRunning(self):
    """
    Return an information concerning the Hub running status.

    @return: True if the hub is running, False otherwise
    @rtype: boolean
    """
    return self._is_running

  def _serve_forever(self):

    if self._web_profile:
      self._web_profile_popup_dialogue = WebProfilePopupDialogue(self._web_profile_requests_result)

    while self._is_running:
      
      try:
        r = w = e = None
        r, w, e = select.select([self._server.socket], [], [], 0.1)
      except:
        pass
      if r:
        self._server.handle_request()

      if self._web_profile:

        try:
          request = self._web_profile_requests_queue.get_nowait()
          self._web_profile_popup_dialogue.showPopup(request)
        except Queue.Empty:
          pass

        try:
          r = w = e = None
          r, w, e = select.select([self._web_profile_server.socket], [], [], 0.01)
          self._web_profile_popup_dialogue.update()
        except:
          pass
        if r:
          self._web_profile_server.handle_request()


  def _notifyShutdown(self):
    msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.shutdown")		
    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        for key in self._mtype2ids[mtype]:
          self._notify_(self._hub_private_key, self._private_keys[key][0],
                        {"samp.mtype":"samp.hub.event.shutdown",
                         "samp.params": {}})

  def _notifyRegister(self, private_key):
    msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.register")		
    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        public_id = self._private_keys[private_key][0]
        for key in self._mtype2ids[mtype]:
          #if key != private_key:
          self._notify(self._hub_private_key, self._private_keys[key][0],
                       {"samp.mtype":"samp.hub.event.register",
                        "samp.params": {"id": public_id}})

  def _notifyUnregister(self, private_key):
    msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.unregister")		
    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        public_id = self._private_keys[private_key][0]
        for key in self._mtype2ids[mtype]:
          if key != private_key:
            self._notify(self._hub_private_key, self._private_keys[key][0],
                         {"samp.mtype":"samp.hub.event.unregister",
                          "samp.params": {"id": public_id}})

  def _notifyMetadata(self, private_key):
    msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.metadata")		
    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        public_id = self._private_keys[private_key][0]
        for key in self._mtype2ids[mtype]:
          #if key != private_key:
          self._notify(self._hub_private_key, self._private_keys[key][0],
                       {"samp.mtype":"samp.hub.event.metadata",
                        "samp.params": {"id": public_id,
                                        "metadata": self._metadata[private_key]}
                        })

  def _notifySubscriptions(self, private_key):
    msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.event.subscriptions")		
    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        public_id = self._private_keys[private_key][0]
        for key in self._mtype2ids[mtype]:
          #if key != private_key:
          self._notify(self._hub_private_key, self._private_keys[key][0],
                       {"samp.mtype":"samp.hub.event.subscriptions",
                        "samp.params": {"id": public_id,
                                        "subscriptions": self._id2mtypes[private_key]}
                        })

  def _notifyDisconnection(self, private_key):

    def _xmlrpc_call_disconnect(endpoint, private_key, hub_public_id, message):
      try:
        endpoint.samp.client.receiveNotification(private_key, hub_public_id, message)
      except:
        pass

    msubs = SAMPHubServer.getMTypeSubtypes("samp.hub.disconnect")
    public_id = self._private_keys[private_key][0]
    endpoint = self._xmlrpcEndpoints[public_id][1]

    for mtype in msubs:
      if self._mtype2ids.has_key(mtype) and private_key in self._mtype2ids[mtype]:
        try:
          self._log.debug("notify disconnection to %s" % (public_id))
          threading.Thread(target = _xmlrpc_call_disconnect,
                           args = (endpoint, private_key, self._hub_public_id,
                                   {"samp.mtype":"samp.hub.disconnect",
                                    "samp.params": {"reason": "Timeout expired!"}})).start()
        except:
          self._log.warning("disconnection notification to client %s failed\n" % (public_id))


  def _ping(self):
    self._updateLastActivityTime()
    self._log.debug("ping")
    return "1"
  
  def _query_by_metadata(self, key, value):
    public_id_list = []
    for private_id in self._metadata:
      if key in self._metadata[private_id]:
        if self._metadata[private_id][key] == value:
          public_id_list.append(self._private_keys[private_id][0])
          
    return public_id_list

  def _setXmlrpcCallback(self, private_key, xmlrpc_addr):
    self._updateLastActivityTime(private_key)
    if self._private_keys.has_key(private_key):
      if private_key == self._hub_private_key:
        self._xmlrpcEndpoints[self._private_keys[private_key][0]] = (xmlrpc_addr, _HubAsClient(self._hubAsClientRequestHandler))
        return ""
      # Dictionary stored with the public id
      self._log.debug("setXmlrpcCallback: %s %s" % (private_key, xmlrpc_addr))
      server_proxy_pool = None
      if SSL_SUPPORT and xmlrpc_addr[0:5] == "https":
        server_proxy_pool = ServerProxyPool(self._pool_size, xmlrpclib.ServerProxy,
                                            xmlrpc_addr, transport = SafeTransport(key_file = self._keyfile,
                                                                                   cert_file = self._certfile,
                                                                                   cert_reqs = self._cert_reqs,
                                                                                   ca_certs = self._ca_certs,
                                                                                   ssl_version = ssl.PROTOCOL_SSLv3),
                                            allow_none=1)
      else:
        server_proxy_pool = ServerProxyPool(self._pool_size, xmlrpclib.ServerProxy, \
                                            xmlrpc_addr, allow_none=1)

      self._xmlrpcEndpoints[self._private_keys[private_key][0]] = (xmlrpc_addr, server_proxy_pool)
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    return ""


  def _perform_standard_register(self):

    self._thread_lock.acquire()
    private_key, public_id = self._getNewIds()
    self._thread_lock.release()
    self._private_keys[private_key] = (public_id, time.time())
    self._updateLastActivityTime(private_key)
    self._notifyRegister(private_key)
    self._log.debug("register: private-key = %s and self-id = %s" % (private_key, public_id))
    return {"samp.self-id": public_id, \
            "samp.private-key": private_key, \
            "samp.hub-id": self._hub_public_id}


  def _register(self, secret):
    self._updateLastActivityTime()
    if secret == self._hub_secret:
      return self._perform_standard_register()
    else:
      #return {"samp.self-id": "", "samp.private-key": "", "samp.hub-id": ""}
      raise SAMPProxyError(7, "Bad secret code")

  def _getNewIds(self):

    now = datetime.datetime.utcnow().isoformat()
    host = self._host_name
    user = getpass.getuser()
    private_key = hashlib.md5(str(random.randint(0, 99999999999999999999)) \
                              + self._hub_secret + host + user + now).hexdigest()
    self._client_id_counter += 1
    public_id = 'cli#hub'
    if self._client_id_counter > 0:
      public_id = "cli#%d" % (self._client_id_counter)

    return private_key, public_id


  def _unregister(self, private_key):

    self._updateLastActivityTime()

    public_key = ""

    self._notifyUnregister(private_key)

    self._thread_lock.acquire()

    if self._private_keys.has_key(private_key):
      public_key = self._private_keys[private_key][0]
      del self._private_keys[private_key]
    else:
      self._thread_lock.release()
      return ""

    if self._metadata.has_key(private_key):
      del self._metadata[private_key]

    if self._id2mtypes.has_key(private_key):
      del self._id2mtypes[private_key]

    for mtype in self._mtype2ids.keys():
      if private_key in self._mtype2ids[mtype]:
        self._mtype2ids[mtype].remove(private_key)

    if self._xmlrpcEndpoints.has_key(public_key):
      del self._xmlrpcEndpoints[public_key]

    if self._client_activity_time.has_key(private_key):
      del self._client_activity_time[private_key]

    if self._web_profile:
      if self._web_profile_callbacks.has_key(private_key):
        del self._web_profile_callbacks[private_key]
      self._web_profile_server.remove_client(private_key)

    self._thread_lock.release()

    self._log.debug("unregister %s (%s)" % (public_key, private_key))

    return ""

  def _declareMetadata(self, private_key, metadata):
    self._updateLastActivityTime(private_key)
    if self._private_keys.has_key(private_key):
      self._log.debug("declareMetadata: private-key = %s metadata = %s" % (private_key, str(metadata)))
      self._metadata[private_key] = metadata
      self._notifyMetadata(private_key)
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)
    return ""

  def _getMetadata(self, private_key, client_id):
    self._updateLastActivityTime(private_key)
    if self._private_keys.has_key(private_key):
      client_private_key = self._getPrivateKeyFromPublicId(client_id)
      self._log.debug("getMetadata: private-key = %s client-id = %s" %  \
                      (private_key, client_id))
      if client_private_key != None:
        if self._metadata.has_key(client_private_key):
          self._log.debug("--> metadata = %s" % self._metadata[client_private_key])
          return self._metadata[client_private_key]
        else:
          return {}
      else:
        self._log.warning("Invalid client ID %s" % client_id)
        raise SAMPProxyError(6, "Invalid client ID")
    else:
      self._log.warning("Private-key %s expired or invalid." % private_key)
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _declareSubscriptions(self, private_key, mtypes):

    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):

      self._log.debug("declareSubscriptions: private-key = %s mtypes = %s" % (private_key, str(mtypes)))

      # remove subscription to previous mtypes
      if self._id2mtypes.has_key(private_key):

        prev_mtypes = self._id2mtypes[private_key]

        for mtype in prev_mtypes:
          try:
            self._mtype2ids[mtype].remove(private_key)
          except:
            pass

      self._id2mtypes[private_key] = copy.deepcopy(mtypes)

      # remove duplicated MType for wildcard overwriting
      original_mtypes = copy.deepcopy(mtypes)

      for mtype in original_mtypes:
        if mtype.endswith("*"):
          for mtype2 in original_mtypes:
            if mtype2.startswith(mtype[:-1]) and \
               mtype2 != mtype:
              if mtypes.has_key(mtype2): del(mtypes[mtype2])

      self._log.debug("declareSubscriptions: subscriptions accepted from %s => %s" % (private_key, str(mtypes)))

      for mtype in mtypes:

        if self._mtype2ids.has_key(mtype):
          if not private_key in self._mtype2ids[mtype]:
            self._mtype2ids[mtype].append(private_key)
        else:
          self._mtype2ids[mtype] = [private_key]

      self._notifySubscriptions(private_key)

    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    return ""


  def _getSubscriptions(self, private_key, client_id):

    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      client_private_key = self._getPrivateKeyFromPublicId(client_id)
      if client_private_key != None:
        if self._id2mtypes.has_key(client_private_key):
          self._log.debug("getSubscriptions: client-id = %s mtypes = %s" %\
                          (client_id, str(self._id2mtypes[client_private_key])))
          return self._id2mtypes[client_private_key]
        else:
          self._log.debug("getSubscriptions: client-id = %s mtypes = missing" % client_id)
          return {}
      else:
        raise SAMPProxyError(6, "Invalid client ID")
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)


  def _getRegisteredClients(self, private_key):

    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      reg_clients = []
      for pkey in self._private_keys.keys():
        if pkey != private_key:
          reg_clients.append(self._private_keys[pkey][0])
      self._log.debug("getRegisteredClients: private_key = %s clients = %s" % (private_key, reg_clients))
      return reg_clients
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _getSubscribedClients(self, private_key, mtype):

    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      sub_clients = {}

      for pkey in self._private_keys.keys():
        if pkey != private_key and self._isSubscribed(pkey, mtype):
          sub_clients[self._private_keys[pkey][0]] =  {}

      self._log.debug("getSubscribedClients: private_key = %s mtype = %s clients = %s" % \
                      (private_key, mtype, sub_clients))
      return sub_clients
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)


  @staticmethod
  def getMTypeSubtypes(mtype):
    """
    Return a list containing all the possible wildcarded subtypes of MType. Example:

    >>> import sampy
    >>> sampy.SAMPHubServer.getMTypeSubtypes("samp.app.ping")
    ['samp.app.ping', 'samp.app.*', 'samp.*', '*']

    @param mtype: the MType to be parsed
    @type mtype: string
    @return: a list of subtypes
    @rtype: list
    """

    subtypes = []

    msubs = mtype.split(".")
    indexes = range(len(msubs))
    indexes.reverse()
    indexes.append(-1)

    for i in indexes:
      tmp_mtype = string.join(msubs[:i+1], ".")
      if tmp_mtype != mtype:
        if tmp_mtype != "":
          tmp_mtype = tmp_mtype + ".*"
        else:
          tmp_mtype = "*"
      subtypes.append(tmp_mtype)

    return 	subtypes

  def _isSubscribed(self, private_key, mtype):

    subscribed = False

    msubs = SAMPHubServer.getMTypeSubtypes(mtype)

    for msub in msubs:
      if self._mtype2ids.has_key(msub):
        if private_key in self._mtype2ids[msub]:
          subscribed = True

    return subscribed


  def _notify(self, private_key, recipient_id, message):
    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      if self._isSubscribed(self._getPrivateKeyFromPublicId(recipient_id),
                            message["samp.mtype"]) == False:
        raise SAMPProxyError(2, "Client %s not subscribed to MType %s" % (recipient_id, message["samp.mtype"]))

      threading.Thread(target = self._notify_, args = (private_key, recipient_id, message)).start()
      return {}
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _notify_(self, sender_private_key, recipient_public_id, message):
    if self._private_keys.has_key(sender_private_key):
      sender_public_id = self._private_keys[sender_private_key][0]
      try:
        self._log.debug("notify %s from %s to %s" % (message["samp.mtype"], sender_public_id, recipient_public_id))
        recipient_private_key = self._getPrivateKeyFromPublicId(recipient_public_id)
        if recipient_private_key != None:
          for attempt in xrange(10):
            if self._is_running:
              try:
                if self._web_profile and recipient_private_key in self._web_profile_callbacks:
                  # Web Profile
                  self._web_profile_callbacks[recipient_private_key].put({"samp.methodName": "receiveNotification",
                                                                          "samp.params": [sender_public_id, message]})
                  return
                # Standard Profile
                self._xmlrpcEndpoints[recipient_public_id][1].samp.client.receiveNotification(\
                  recipient_private_key, sender_public_id, message)
                break
              except:
                err = StringIO.StringIO()
                traceback.print_exc(file=err)
                txt = err.getvalue()
                self._log.debug("%s XML-RPC endpoint error (attempt %d): \n%s" % (recipient_public_id, attempt + 1, txt))
                time.sleep(0.01)
        else:
          raise SAMPProxyError(6, "Invalid client ID")
      except:
        self._log.warning("%s notification from client %s to client %s failed\n" % (message["samp.mtype"], sender_public_id, recipient_public_id))

  def _notifyAll(self, private_key, message):
    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      if not message.has_key("samp.mtype"):
        raise SAMPProxyError(3, "samp.mtype keyword is missing")
      recipient_ids = self._notifyAll(private_key, message)
      return recipient_ids
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _notifyAll(self, sender_private_key, message):

    recipient_ids = []
    msubs = SAMPHubServer.getMTypeSubtypes(message["samp.mtype"])

    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        for key in self._mtype2ids[mtype]:
          if key != sender_private_key:
            _recipient_id = self._private_keys[key][0]
            recipient_ids.append(_recipient_id)
            threading.Thread(target=self._notify, 
                             args=(sender_private_key,
                                   _recipient_id, message)
                             ).start()

    return recipient_ids

  def _call(self, private_key, recipient_id, msg_tag, message):
    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      if self._isSubscribed(self._getPrivateKeyFromPublicId(recipient_id),
                            message["samp.mtype"]) == False:
        raise SAMPProxyError(2, "Client %s not subscribed to MType %s" % (recipient_id, message["samp.mtype"]))
      public_id = self._private_keys[private_key][0]
      msg_id = self._getNewHubMsgId(public_id, msg_tag)
      threading.Thread(target = self._call_, args = (private_key, public_id, recipient_id, msg_id, message)).start()
      return msg_id
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _call_(self, sender_private_key, sender_public_id, recipient_public_id, msg_id, message):
    if self._private_keys.has_key(sender_private_key):
      try:
        self._log.debug("call %s from %s to %s (%s)" % (msg_id.split(";;")[0], sender_public_id, recipient_public_id, message["samp.mtype"]))
        recipient_private_key = self._getPrivateKeyFromPublicId(recipient_public_id)
        if recipient_private_key != None:
          for attempt in xrange(10):
            if self._is_running:
              try:
                if self._web_profile and recipient_private_key in self._web_profile_callbacks:
                  # Web Profile
                  self._web_profile_callbacks[recipient_private_key].put({"samp.methodName": "receiveCall",
                                                                          "samp.params": [sender_public_id, msg_id, message]})
                  return
                # Standard Profile
                self._xmlrpcEndpoints[recipient_public_id][1].samp.client.receiveCall(\
                  recipient_private_key, sender_public_id, msg_id, message)
                break
              except:
                err = StringIO.StringIO()
                traceback.print_exc(file=err)
                txt = err.getvalue()
                self._log.debug("%s XML-RPC endpoint error (attempt %d): \n%s" % (recipient_public_id, attempt + 1, txt))
                time.sleep(0.01)
        else:
          raise SAMPProxyError(6, "Invalid client ID")
      except:
        self._log.warning("%s call %s from client %s to client %s failed\n" % (message["samp.mtype"], msg_id.split(";;")[0], sender_public_id, recipient_public_id))

  def _callAll(self, private_key, msg_tag, message):
    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      if not message.has_key("samp.mtype"):
        raise SAMPProxyError(3, "samp.mtype keyword is missing in message tagged as %s" % msg_tag)

      public_id = self._private_keys[private_key][0]
      msg_id = self._callAll_(private_key, public_id, msg_tag, message)
      return msg_id
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _callAll_(self, sender_private_key, sender_public_id, msg_tag, message):

    msg_id = {}
    msubs = SAMPHubServer.getMTypeSubtypes(message["samp.mtype"])

    for mtype in msubs:
      if self._mtype2ids.has_key(mtype):
        for key in self._mtype2ids[mtype]:
          if key != sender_private_key:
            _msg_id = self._getNewHubMsgId(sender_public_id, msg_tag)
            receiver_public_id = self._private_keys[key][0]
            msg_id[receiver_public_id] = _msg_id
            threading.Thread(target=self._call_, 
                             args=(sender_private_key, sender_public_id,
                                   receiver_public_id, _msg_id, message)
                             ).start()
    return msg_id

  def _callAndWait(self, private_key, recipient_id, message, timeout):
    self._updateLastActivityTime(private_key)

    if self._private_keys.has_key(private_key):
      timeout = int(timeout)

      now = time.time()
      response = {}

      msg_id = self._call(private_key, recipient_id, "sampy::sync::call", message)
      self._sync_msg_ids_heap[msg_id] = None

      while self._is_running:
        if timeout > 0 and time.time() - now >= timeout:
          del(self._sync_msg_ids_heap[msg_id])
          raise SAMPProxyError(1, "Timeout expired!")

        if self._sync_msg_ids_heap[msg_id] != None:
          response = copy.deepcopy(self._sync_msg_ids_heap[msg_id])
          del(self._sync_msg_ids_heap[msg_id])
          break
        time.sleep(0.01)

      return response
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

  def _reply(self, private_key, msg_id, response):
    self._updateLastActivityTime(private_key)
    if self._private_keys.has_key(private_key):
      threading.Thread(target = self._reply_, args = (private_key, msg_id, response)).start()
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)

    return {}

  def _reply_(self, responder_private_key, msg_id, response):

    if self._private_keys.has_key(responder_private_key) and msg_id:
      responder_public_id = self._private_keys[responder_private_key][0]
      counter, hub_public_id, recipient_public_id, \
             recipient_msg_tag = msg_id.split(";;", 3)

      try:
        self._log.debug("reply %s from %s to %s" % (counter, responder_public_id, recipient_public_id))
        if recipient_msg_tag == "sampy::sync::call":
          if msg_id in self._sync_msg_ids_heap.keys():
            self._sync_msg_ids_heap[msg_id] = response
        else:
          recipient_private_key = self._getPrivateKeyFromPublicId(recipient_public_id)
          if recipient_private_key != None:
            for attempt in xrange(10):
              if self._is_running:
                try:
                  if self._web_profile and recipient_private_key in self._web_profile_callbacks:
                    # Web Profile
                    self._web_profile_callbacks[recipient_private_key].put({"samp.methodName": "receiveResponse",
                                                                            "samp.params": [responder_public_id, \
                                                                                            recipient_msg_tag, response]})
                    return
                  # Standard Profile
                  self._xmlrpcEndpoints[recipient_public_id][1].samp.client.receiveResponse(\
                    recipient_private_key, responder_public_id, recipient_msg_tag, response)
                  break
                except:
                  err = StringIO.StringIO()
                  traceback.print_exc(file=err)
                  txt = err.getvalue()
                  self._log.debug("%s XML-RPC endpoint error (attempt %d): \n%s" % (recipient_public_id, attempt + 1, txt))
                  time.sleep(0.01)
          else:
            raise SAMPProxyError(6, "Invalid client ID")
      except:
        self._log.warning("%s reply from client %s to client %s failed\n" % (recipient_msg_tag, responder_public_id, recipient_public_id))

  def _getPrivateKeyFromPublicId(self, public_id):

    for private_key in self._private_keys.keys():
      if self._private_keys[private_key][0] == public_id:
        return private_key
    return None

  def _getNewHubMsgId(self, sender_public_id, sender_msg_id):
    self._thread_lock.acquire()
    self._hub_msg_id_counter += 1
    self._thread_lock.release()
    return "msg#%d;;%s;;%s;;%s" % \
           (self._hub_msg_id_counter, self._hub_public_id, \
            sender_public_id, sender_msg_id)

  def _updateLastActivityTime(self, private_key = None):
    self._thread_lock.acquire()
    self._last_activity_time = time.time()
    if private_key != None:
      self._client_activity_time[private_key] = time.time()
    self._thread_lock.release()

  def _receiveNotification(self, private_key, sender_id, message):
    return ""

  def _receiveCall(self, private_key, sender_id, msg_id, message):
    if private_key == self._hub_private_key:
      
      if message.has_key("samp.mtype") and message["samp.mtype"] == "samp.app.ping":
        self._reply(self._hub_private_key, msg_id, {"samp.status": SAMP_STATUS_OK, "samp.result": {}})
        
      elif message.has_key("samp.mtype") and \
         (message["samp.mtype"] == "x-samp.query.by-meta" or \
          message["samp.mtype"] == "samp.query.by-meta"):
        
        ids_list = self._query_by_metadata(message["samp.params"]["key"], message["samp.params"]["value"])
        self._reply(self._hub_private_key, msg_id, 
                    {"samp.status": SAMP_STATUS_OK, 
                     "samp.result": {"ids": ids_list}})
      
      return ""
    else:
      return ""

  def _receiveResponse(self, private_key, responder_id, msg_tag, response):
    return ""


  def _web_profile_register(self, identity_info, client_address = ("uknown", 0), origin = "unknown"):
    
    self._updateLastActivityTime()

    if not client_address[0] in ["localhost", "127.0.0.1"]:
      raise SAMPProxyError(403, "Request of registration rejected by the Hub.")

    if not origin:
      origin = "unknown"

    if isinstance(identity_info, dict):
      # an old version of the protocol provided just a string with the app name
      if "samp.name" not in identity_info.keys():
        raise SAMPProxyError(403, "Request of registration rejected by the Hub (application name not provided).")
    
    # Red semaphore for the other threads
    self._web_profile_requests_semaphore.put("wait")
    # Set the request to be displayed for the current thread
    self._web_profile_requests_queue.put((identity_info, client_address, origin))
    # Get the popup dialogue response
    response = self._web_profile_requests_result.get()
    # OK, semaphore green
    self._web_profile_requests_semaphore.get()

    if response:
      register_map = self._perform_standard_register()
      register_map["samp.url-translator"] = "http://localhost:21012/translator/%s?ref=" % register_map["samp.private-key"]
      self._web_profile_server.add_client(register_map["samp.private-key"])
      return register_map
    else:
      raise SAMPProxyError(403, "Request of registration rejected by the user.")


  def _web_profile_allowReverseCallbacks(self, private_key, allow):
    self._updateLastActivityTime()
    if self._private_keys.has_key(private_key):
      if allow == "0":
        if private_key in self._web_profile_callbacks:
          del self._web_profile_callbacks[private_key]
      else:
        self._web_profile_callbacks[private_key] = Queue.Queue()
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)
    return ""

  def _web_profile_pullCallbacks(self, private_key, timeout_secs):
    self._updateLastActivityTime()
    if self._private_keys.has_key(private_key):
      callback = []
      try:
        while self._is_running:
          item_queued = self._web_profile_callbacks[private_key].get(block=True, timeout=int(timeout_secs))
          callback.append(item_queued)
          if self._web_profile_callbacks[private_key].empty():
            break
      except Queue.Empty:
        pass
      return callback
    else:
      raise SAMPProxyError(5, "Private-key %s expired or invalid." % private_key)



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


class SAMPHubProxy(object):
  """
  Proxy class useful to simplify the client interaction with a SAMP Hub.
  """

  def __init__(self):
    self.proxy = None
    self._connected = False
    self.lockfile = {}

  def isConnected(self):
    """
    Testing method to verify the proxy connection with a running Hub.

    @return: return True if the proxy is connected to a Hub, False otherwise
    @rtype: boolean
    """
    return self._connected


  @staticmethod
  def getRunningHubs():
    """
    Return a dictionary containing the lock-file contents of all the currently
    running hubs (single and/or multiple mode). The dictionary format is:

    C{{<lock-file>: {<token-name>: <token-string>, ...}, ...}}

    where C{<lock-file>} is the lock-file name, C{<token-name>} and C{<token-string>}
    are the lock-file tokens (name and content).

    @return: the lock-file contents of all the currently running hubs.
    @rtype: dictionary

    """

    hubs = {}
    lockfilename = ""

    # HUB SINGLE INSTANCE MODE

    # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
    if os.environ.has_key("SAMP_HUB"):
      # For the time being I assume just the std profile supported.
      if os.environ["SAMP_HUB"].startswith("std-lockurl:"):
        lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:"):]
    else:
      if os.environ.has_key("HOME"):
        # UNIX
        lockfilename = os.path.join(os.environ["HOME"], ".samp")
      else:
        # Windows
        lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")

    hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)

    if hub_is_running:
      hubs[lockfilename] = lockfiledict

    # HUB MULTIPLE INSTANCE MODE

    lockfiledir = ""

    if os.environ.has_key("HOME"):
      # UNIX
      lockfiledir = os.path.join(os.environ["HOME"], ".samp-1")
    else:
      # Windows
      lockfiledir = os.path.join(os.environ["USERPROFILE"], ".samp-1")

    if os.path.isdir(lockfiledir):
      for filename in os.listdir(lockfiledir):
        if re.match('samp\\-hub\\-\d+\\-\d+', filename) != None:
          lockfilename = os.path.join(lockfiledir, filename)
          hub_is_running, lockfiledict = SAMPHubServer.checkRunningHub(lockfilename)
          if hub_is_running:
            hubs[lockfilename] = lockfiledict

    return hubs

  def connect(self, hub_params = None, user = None, password = None,
              key_file=None, cert_file=None, cert_reqs=0,
              ca_certs=None, ssl_version=1, pool_size=20):
    """
    Connect to the current SAMP Hub. If a SAMP Hub is not running or refuses the connection,
    then a L{SAMPHubError} is raised.

    @param hub_params: Optional dictionary containig the lock-file content of the Hub
    with which to connect. This dictionary has the form C{{<token-name>: <token-string>, ...}}.
    @type hub_params: dictionary

    @param user: In case of Basic Authenticated connections, C{user} specifies the user name.
    @type user: string

    @param password: In case of Basic Authenticated connections, C{password} specifies the user
    password.
    @type password: string

    @param key_file: Set the file containing the private key for SSL connections. If the
    certificate file (C{certfile}) contains the private key, then C{keyfile} can be omitted.
    @type key_file: string

    @param cert_file: Specify the file which contains a certificate to be used to identify the
    local side of the secure connection.
    @type cert_file: string

    @param cert_reqs: The parameter C{cert_reqs} specifies whether a certificate is required
    from the server side of the connection, and whether it will be validated if provided. It
    must be one of the three values L{ssl.CERT_NONE} (certificates ignored), L{ssl.CERT_OPTIONAL}
    (not required, but validated if provided), or L{ssl.CERT_REQUIRED} (required and validated).
    If the value of this parameter is not L{ssl.CERT_NONE}, then the C{ca_certs} parameter must
    point	to a file of CA certificates.
    @type cert_reqs: int

    @param ca_certs: The C{ca_certs} file contains a set of concatenated \"Certification Authority\" 
    certificates, which are used to validate the certificate passed from the server end of the 
    connection.
    @type ca_certs: string

    @param ssl_version: The C{ssl_version} option specifies which version of the SSL protocol to use.
    Typically, the server chooses	a particular protocol version, and the client must adapt to the
    server's choice. Most of the versions are	not interoperable with the other versions. If not
    specified the default SSL version is  L{ssl.PROTOCOL_SSLv3}. This version provides the most
    compatibility with other versions server side. Other SSL protocol versions are:                        
    L{ssl.PROTOCOL_SSLv2}, L{ssl.PROTOCOL_SSLv23} and L{ssl.PROTOCOL_TLSv1}.
    @type ssl_version: int

    @param pool_size: The number of socket connections opened to communicate with the Hub.
    @type callable: int
    """

    self._connected = False
    self.lockfile = {}

    if hub_params == None:
      hubs = SAMPHubProxy.getRunningHubs()
      if len(hubs.keys()) > 0:
        # Use Single instance hub by default
        lockfilename = ""
        # CHECK FOR SAMP_HUB ENVIRONMENT VARIABLE
        if os.environ.has_key("SAMP_HUB"):
          # For the time being I assume just the std profile supported.
          if os.environ["SAMP_HUB"].startswith("std-lockurl:"):
            lockfilename = os.environ["SAMP_HUB"][len("std-lockurl:"):]
          else:
            raise SAMPHubError("SAMP Hub profile not supported.")
        else:
          if os.environ.has_key("HOME"):
            # UNIX
            lockfilename = os.path.join(os.environ["HOME"], ".samp")
          else:
            # Windows
            lockfilename = os.path.join(os.environ["USERPROFILE"], ".samp")
        hub_params = hubs[lockfilename]
      else:
        raise SAMPHubError("Unable to find a running SAMP Hub.")

    try:

      url = hub_params["samp.hub.xmlrpc.url"].replace("\\", "")

      # URL formatting for Basic Authentication parameters
      if user != None and password != None:
        trans, addr = url.split("://")
        url = "%s://%s:%s@%s" % (trans, user, password, addr)

      if SSL_SUPPORT and url[0:5] == "https":
        self.proxy = ServerProxyPool(pool_size, xmlrpclib.ServerProxy,
                                     url, transport = SafeTransport(key_file, cert_file, cert_reqs,
                                                                    ca_certs, ssl_version),
                                     allow_none=1)
      else:
        self.proxy = ServerProxyPool(pool_size, xmlrpclib.ServerProxy, url, allow_none=1)

      self.proxy.samp.hub.ping()

      self.lockfile = copy.deepcopy(hub_params)
      self._connected = True

    except xmlrpclib.ProtocolError, p:
      # 401 Unauthorized
      if p.errcode == 401:
        raise SAMPHubError("Unauthorized access. Basic Authentication required or failed.")
      else:
        raise SAMPHubError("Protocol Error %d: %s" % (p.errcode, p.errmsg))
    except:
      err = StringIO.StringIO()
      traceback.print_exc(file=err)
      txt = err.getvalue()
      if SSL_SUPPORT:
        if sys.exc_info()[0] == ssl.SSLError:
          raise SAMPHubError("SSL Error: %s" % sys.exc_info()[1])
        else:
          raise SAMPHubError("SAMP Hub connection refused.\n " + txt)
      else:
        raise SAMPHubError("SAMP Hub connection refused.\n" + txt)


  def disconnect(self):
    """
    Disconnect from the current SAMP Hub.
    """

    self.proxy = None
    self._connected = False
    self.lockfile = {}

  def ping(self):
    """
    Proxy to C{ping} SAMP Hub method (Standard Profile only)
    """
    return self.proxy.samp.hub.ping()

  def setXmlrpcCallback(self, private_key, xmlrpc_addr):
    """
    Proxy to C{setXmlrpcCallback} SAMP Hub method (Standard Profile only)
    """
    return self.proxy.samp.hub.setXmlrpcCallback(private_key, xmlrpc_addr)

  def register(self, secret):
    """
    Proxy to C{register} SAMP Hub method
    """
    return self.proxy.samp.hub.register(secret)

  def unregister(self, private_key):
    """
    Proxy to C{unregister} SAMP Hub method
    """
    return self.proxy.samp.hub.unregister(private_key)

  def declareMetadata(self, private_key, metadata):
    """
    Proxy to C{declareMetadata} SAMP Hub method
    """
    return self.proxy.samp.hub.declareMetadata(private_key, metadata)

  def getMetadata(self, private_key, client_id):
    """
    Proxy to C{getMetadata} SAMP Hub method
    """
    return self.proxy.samp.hub.getMetadata(private_key, client_id)

  def declareSubscriptions(self, private_key, subscriptions):
    """
    Proxy to C{declareSubscriptions} SAMP Hub method
    """
    return self.proxy.samp.hub.declareSubscriptions(private_key, subscriptions)

  def getSubscriptions(self, private_key, client_id):
    """
    Proxy to C{getSubscriptions} SAMP Hub method
    """
    return self.proxy.samp.hub.getSubscriptions(private_key, client_id)

  def getRegisteredClients(self, private_key):
    """
    Proxy to C{getRegisteredClients} SAMP Hub method
    """
    return self.proxy.samp.hub.getRegisteredClients(private_key)

  def getSubscribedClients(self, private_key, mtype):
    """
    Proxy to C{getSubscribedClients} SAMP Hub method
    """
    return self.proxy.samp.hub.getSubscribedClients(private_key, mtype)

  def notify(self, private_key, recipient_id, message):
    """
    Proxy to C{notify} SAMP Hub method
    """
    # Add user in Basic Authentication case
    return self.proxy.samp.hub.notify(private_key, recipient_id, message)

  def notifyAll(self, private_key, message):
    """
    Proxy to C{notifyAll} SAMP Hub method
    """
    return self.proxy.samp.hub.notifyAll(private_key, message)

  def call(self, private_key, recipient_id, msg_tag, message):
    """
    Proxy to C{call} SAMP Hub method
    """
    return self.proxy.samp.hub.call(private_key, recipient_id, msg_tag, message)

  def callAll(self, private_key, msg_tag, message):
    """
    Proxy to C{callAll} SAMP Hub method
    """
    return self.proxy.samp.hub.callAll(private_key, msg_tag, message)

  def callAndWait(self, private_key, recipient_id, message, timeout):
    """
    Proxy to C{callAndWait} SAMP Hub method. If timeout expires a 
    L{SAMPProxyError} instance is raised.
    """
    return self.proxy.samp.hub.callAndWait(private_key, recipient_id, message, timeout)

  def reply(self, private_key, msg_id, response):
    """
    Proxy to C{reply} SAMP Hub method
    """
    return self.proxy.samp.hub.reply(private_key, msg_id, response)


_THREAD_STARTED_COUNT = 0

class SAMPClient(object):
  """
  Utility class which provides facilities to create and manage a SAMP compliant
  XML-RPC server that acts as SAMP callable client application. 
  """

  def __init__(self, hub, name = None, description=None, metadata = None, \
               addr = None, port = 0, https = False, key_file=None, cert_file=None, \
               cert_reqs=0, ca_certs=None, ssl_version=2, callable = True):
    """
    L{SAMPClient} constructor.

    @param hub: an instance of L{SAMPHubProxy} to be used for messaging with the SAMP Hub.
    @type hub: L{SAMPHubProxy}
    
    @param name: (optional) a string containing the client name
    (corresponding to C{samp.name} metadata keyword).
    @type name: string
    
    @param description: (optional) a string containing the client description
    (corresponding to C{samp.description.text} metadata keyword).
    @type description: string

    @param metadata: (optional) a dictionary containing the client 
    application metadata in the standard SAMP format
    @type metadata: dict

    @param addr: Listening address (or IP)
    @type addr: string

    @param port: (optional) the listening XML-RPC server socket port
    @type port: int

    @param https: set the callable client running on a Secure Sockets Layer connection (HTTPS).
    By default SSL is desabled.
    @type https: boolean

    @param key_file: Set the file containing the private key for SSL connections. If the
    certificate file (C{cert_file}) contains the private key, then C{key_file} can be omitted.
    @type key_file: string

    @param cert_file: Specify the file which contains a certificate to be used to identify the
    local side of the secure connection.
    @type cert_file: string

    @param cert_reqs: The parameter C{cert_reqs} specifies whether a certificate is required
    from the Hub side of the connection, and whether it will be validated if provided. It
    must be one of the three values L{ssl.CERT_NONE} (certificates ignored), L{ssl.CERT_OPTIONAL}
    (not required, but validated if provided), or L{ssl.CERT_REQUIRED} (required and validated).
    If the value of this parameter is not L{ssl.CERT_NONE}, then the C{ca_certs} parameter must
    point	to a file of CA certificates.
    @type cert_reqs: int

    @param ca_certs: The C{ca_certs} file contains a set of concatenated \"Certification Authority\" 
    certificates, which are used to validate the certificate passed from the Hub end of the 
    connection.
    @type ca_certs: string

    @param ssl_version: The C{ssl_version} option specifies which version of the SSL protocol to use.
    Typically, the server chooses	a particular protocol version, and the client must adapt to the
    server's choice. Most of the versions are	not interoperable with the other versions. If not
    specified the default SSL version is  L{ssl.PROTOCOL_SSLv23}. This version provides the most
    compatibility with other versions Hub side. Other SSL protocol versions are:                        
    L{ssl.PROTOCOL_SSLv2}, L{ssl.PROTOCOL_SSLv3} and L{ssl.PROTOCOL_TLSv1}.
    @type ssl_version: int

    @param callable: Specify whether the client is a callable client or not
    @type callable: boolean

    """

    # GENERAL
    self._thread = None
    self._is_running = False
    
    if metadata == None: metadata = {}
    
    if name != None:
      metadata["samp.name"] = name
      
    if description != None:
      metadata["samp.description.text"] = description
    
    self._metadata = metadata
    
    self._addr = addr
    self._port = port
    self._xmlrpcAddr = None
    self._callable = callable

    # HUB INTERACTION
    self.client = None
    self._public_id = None
    self._private_key = None
    self._hub_id = None
    self._notification_bindings = {}
    self._call_bindings = {"samp.app.ping": [self._ping, {}],
                           "client.env.get": [self._client_env_get, {}]}
    self._response_bindings = {}
    

    self._host_name = "127.0.0.1"
    try:
      if internet_on():
        self._host_name = socket.getfqdn()
    except:
      pass

    self.hub = hub

    if self._callable:

      if SSL_SUPPORT and https:
        self.client = SecureXMLRPCServer((self._addr or self._host_name, self._port), 
                                         keyfile, certfile, cert_reqs, ca_certs, ssl_version,
                                         self._log, logRequests = False, allow_none = True)
      else:
        self.client = ThreadingXMLRPCServer((self._addr or self._host_name,
                                             self._port), logRequests = False, allow_none = True)

      self.client.register_introspection_functions()
      self.client.register_function(self.receiveNotification, 'samp.client.receiveNotification')
      self.client.register_function(self.receiveCall, 'samp.client.receiveCall')
      self.client.register_function(self.receiveResponse, 'samp.client.receiveResponse')

      if self._port == 0:
        self._port = self.client.socket.getsockname()[1]

      self._xmlrpcAddr = "http://%s:%s" % (self._addr or \
                                           self._host_name, \
                                           self._port)

  def __del__(self):
    self.stop()

  def start(self):
    """
    Start the client in a non-blocking way.
    """
    global _THREAD_STARTED_COUNT
    _THREAD_STARTED_COUNT += 1
    self._is_running = True
    self._run_client()

  def stop(self, timeout=0.1):
    """
    Stop the client.

    @param timeout: timeout after wich the client terminates even if the threading is still alive.
    @type timeout: float
    """

    self._is_running = False
    if self._thread is not None:
      self._thread.join(timeout)
      self._thread = None


  def isRunning(self):
    """
    Return an information concerning the client running status.

    @return: True if the client is running, False otherwise
    @rtype: boolean
    """
    return self._is_running != None


  def _run_client(self):
    if self._callable:
      self._thread = threading.Thread(target = self._serve_forever)
      self._thread.setDaemon(True)
      self._thread.start()

  def _serve_forever(self):
    while self._is_running:
      try:
        r = w = e = None
        r, w, e = select.select([self.client.socket], [], [], 0.1)
      except:
        pass
      if r:
        self.client.handle_request()

  def _ping(self, private_key, sender_id, msg_id, msg_mtype, msg_params, message):
    self.hub.reply(private_key, msg_id, {"samp.status": SAMP_STATUS_OK, "samp.result": {}})

  def _client_env_get(self, private_key, sender_id, msg_id, msg_mtype, msg_params, message):
    if msg_params["name"] in os.environ:
      self.hub.reply(private_key, msg_id, {"samp.status": SAMP_STATUS_OK,
                                           "samp.result": {"value": os.environ[msg_params["name"]]}})
    else:
      self.hub.reply(private_key, msg_id, {"samp.status": SAMP_STATUS_WARNING,
                                           "samp.result": {"value": ""},
                                           "samp.error": {"samp.errortxt": 
                                                          "Environment variable not defined."}})

  def _handle_notification(self, private_key, sender_id, message):

    if private_key == self.getPrivateKey() and message.has_key("samp.mtype"):

      msg_mtype = message["samp.mtype"]
      del message["samp.mtype"]
      msg_params = message["samp.params"]
      del message["samp.params"]

      msubs = SAMPHubServer.getMTypeSubtypes(msg_mtype)
      for mtype in msubs:
        if self._notification_bindings.has_key(mtype):
          bound_func = self._notification_bindings[mtype][0]
          if (inspect.ismethod(bound_func) and bound_func.im_func.func_code.co_argcount == 6) or \
             (inspect.isfunction(bound_func) and bound_func.func_code.co_argcount == 5):
            bound_func(private_key, sender_id, msg_mtype, msg_params, message)
          else:
            bound_func(private_key, sender_id, None, msg_mtype, msg_params, message)

    return ""

  def receiveNotification(self, private_key, sender_id, message):
    """
    Standard callable client C{receiveNotification} method. This method is
    automatically handled when L{bindReceiveNotification} method is used to bind
    distinct operations to MTypes. In case of a customized callable client
    implementation that inherits from L{SAMPClient} class this method should
    be overwritten. ATTENTION: when overwritten, this method must always return
    a string result (even empty).

    @param private_key: the client private key.
    @type private_key: str

    @param sender_id: the sender public ID.
    @type sender_id: str

    @param message: the message received.
    @type message: dict

    @return: any confirmation string.
    @rtype: str
    """
    return self._handle_notification(private_key, sender_id, message)

  def _handle_call(self, private_key, sender_id, msg_id, message):
    if private_key == self.getPrivateKey() and message.has_key("samp.mtype"):

      msg_mtype = message["samp.mtype"]
      del message["samp.mtype"]
      msg_params = message["samp.params"]
      del message["samp.params"]

      msubs = SAMPHubServer.getMTypeSubtypes(msg_mtype)

      for mtype in msubs:
        if self._call_bindings.has_key(mtype):
          self._call_bindings[mtype][0](private_key, sender_id, msg_id, msg_mtype, msg_params, message)

    return ""

  def receiveCall(self, private_key, sender_id, msg_id, message):
    """
    Standard callable client C{receiveCall} method. This method is
    automatically handled when L{bindReceiveCall} method is used to bind
    distinct operations to MTypes. In case of a customized callable client
    implementation that inherits from L{SAMPClient} class this method should
    be overwritten. ATTENTION: when overwritten, this method must always return
    a string result (even empty).

    @param private_key: the client private key.
    @type private_key: str

    @param sender_id: the sender public ID.
    @type sender_id: str

    @param msg_id: the message ID received.
    @type msg_id: str

    @param message: the message received.
    @type message: dict

    @return: any confirmation string.
    @rtype: str
    """
    return self._handle_call(private_key, sender_id, msg_id, message)

  def _handle_response(self, private_key, responder_id, msg_tag, response):
    if private_key == self.getPrivateKey() and self._response_bindings.has_key(msg_tag):
      self._response_bindings[msg_tag](private_key, responder_id, msg_tag, response)
    return ""

  def receiveResponse(self, private_key, responder_id, msg_tag, response):
    """
    Standard callable client C{receiveResponse} method. This method is
    automatically handled when L{bindReceiveResponse} method is used to bind
    distinct operations to MTypes. In case of a customized callable client
    implementation that inherits from L{SAMPClient} class this method should
    be overwritten. ATTENTION: when overwritten, this method must always return
    a string result (even empty).

    @param private_key: the client private key.
    @type private_key: str

    @param responder_id: the responder public ID.
    @type responder_id: str

    @param msg_tag: the response message tag.
    @type msg_tag: str

    @param response: the response received.
    @type response: dict

    @return: any confirmation string.
    @rtype: str
    """
    return self._handle_response(private_key, responder_id, msg_tag, response)




  def bindReceiveMessage(self, mtype, function, declare = True, metadata = None):
    """Bind a specific MType to a function or class method, being intended for
    a call or a notification.

    The function must be of the form:
    C{def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)}

    where C{private_key} is the client private-key, C{sender_id} argument is the
    notification sender ID, C{msg_id} is the Hub message-id (calls only, otherwise is None),
    C{mtype} is the message MType, C{params} is the message parameter set (content of
    "samp.params") and C{extra} is a dictionary containing any extra message map entry.
    The client is automatically declared subscribed to the MType by default.

    @param mtype: the MType to be catched.
    @type mtype: str

    @param function: the application function to be used when C{mtype} is received.
    @type function: function or class method

    @param declare: specify whether the client must be automatically declared as
    subscribed to the MType (see also L{declareSubscriptions}).
    @type declare: boolean

    @param metadata: an optional map containing additional metadata to declare associated
    with the MType subscribed to (see also L{declareSubscriptions}).
    @type metadata: dict
    """

    self.bindReceiveCall(mtype, function, declare = True, metadata = None)
    self.bindReceiveNotification(mtype, function, declare = True, metadata = None)


  def bindReceiveNotification(self, mtype, function, declare = True, metadata = None):
    """
    Bind a specific MType notification to a function or class method.
    The function must be of the form:

    C{def my_function_or_method(<self,> private_key, sender_id, mtype, params, extra)}

    where C{private_key} is the client private-key, C{sender_id} argument is 
    the notification sender ID, C{mtype} is the message MType, C{params} is 
    the notified message parameter set (content of "samp.params") and C{extra} is a 
    dictionary containing any extra message map entry. The client is
    automatically declared subscribed to the MType by default.

    @param mtype: the MType to be catched.
    @type mtype: str

    @param function: the application function to be used when C{mtype} is received.
    @type function: function or class method

    @param declare: specify whether the client must be automatically declared as
    subscribed to the MType (see alse L{declareSubscriptions}).
    @type declare: boolean

    @param metadata: an optional map containing additional metadata to declare associated
    with the MType subscribed to (see also L{declareSubscriptions}).
    @type metadata: dict
    """
    if self._callable:
      if not metadata:
        metadata = {}
      self._notification_bindings[mtype] = [function, metadata]
      if declare: self._declareSubscriptions()
    else:
      raise SAMPClientError("Client not callable.")

  def bindReceiveCall(self, mtype, function, declare = True, metadata = None):
    """
    Bind a specific MType call to a function or class method.
    The function must be of the form:

    C{def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)}

    where C{private_key} is the client private-key, C{sender_id} argument is the
    notification sender ID, C{msg_id} is the Hub message-id, C{mtype} is the message MType, 
    C{params} is the message parameter set (content of "samp.params") and C{extra} is a 
    dictionary containing any extra message map entry. The client is
    automatically declared subscribed to the MType by default.

    @param mtype: the MType to be catched.
    @type mtype: str

    @param function: the application function to be used when C{mtype} is received.
    @type function: function or class method

    @param declare: specify whether the client must be automatically declared as
    subscribed to the MType (see also L{declareSubscriptions}).
    @type declare: boolean

    @param metadata: an optional map containing additional metadata to declare associated
    with the MType subscribed to (see also L{declareSubscriptions}).
    @type metadata: dict
    """
    if self._callable:
      if not metadata:
        metadata = {}
      self._call_bindings[mtype] = [function, metadata]
      if declare: self._declareSubscriptions()
    else:
      raise SAMPClientError("Client not callable.")

  def bindReceiveResponse(self, msg_tag, function):
    """
    Bind a specific msg-tag response to a function or class method.
    The function must be of the form:

    C{def my_function_or_method(<self,> private_key, responder_id, msg_tag, response)}

    where C{private_key} is the client private-key, C{responder_id} argument is the message
    responder ID, C{msg_tag} is the message-tag provided at call time and C{response} is the
    response received.

    @param msg_tag: the message-tag to be catched.
    @type msg_tag: str

    @param function: the application function to be used when C{msg_tag} is received.
    @type function: function or class method
    """
    if self._callable:
      self._response_bindings[msg_tag] = function
    else:
      raise SAMPClientError("Client not callable.")

  def unbindReceiveNotification(self, mtype, declare = True):
    """
    Remove from the notifications binding table the specified MType and unsubscribe
    the client from it (if required).

    @param mtype: the MType to be removed
    @type mtype: str

    @param declare: specify whether the client must be automatically declared as
    unsubscribed from the MType (see alse L{declareSubscriptions}).
    @type declare: boolean

    """
    if self._callable:
      del self._notification_bindings[mtype]
      if declare: self._declareSubscriptions()
    else:
      raise SAMPClientError("Client not callable.")

  def unbindReceiveCall(self, mtype, declare = True):
    """
    Remove from the calls binding table the specified MType and unsubscribe
    the client from it (if required).

    @param mtype: the MType to be removed
    @type mtype: str

    @param declare: specify whether the client must be automatically declared as
    unsubscribed from the MType (see alse L{declareSubscriptions}).
    @type declare: boolean
    """
    if self._callable:
      del self._call_bindings[mtype]
      if declare: self._declareSubscriptions()
    else:
      raise SAMPClientError("Client not callable.")

  def unbindReceiveResponse(self, msg_tag):
    """
    Remove from the responses binding table the specified message-tag.

    @param msg_tag: the message-tag to be removed
    @type msg_tag: str
    """
    if self._callable:
      del self._response_bindings[msg_tag]
    else:
      raise SAMPClientError("Client not callable.")

  def declareSubscriptions(self, subscriptions = None):
    """
    Declares the MTypes the client wishes to subscribe to, implicitly defined
    with the MType binding methods L{bindReceiveNotification} and L{bindReceiveCall}.
    An optional C{subscriptions} map can be added to the final map passed to 
    the L{SAMPHubProxy.declareSubscriptions} operation.

    @param subscriptions: an optional map containing the list of MTypes to subscribe to,
    with the same format of the C{subscriptions} map passed to the
    L{SAMPHubProxy.declareSubscriptions} operation.
    @type subscriptions: dict
    """
    if self._callable:
      self._declareSubscriptions(subscriptions)
    else:
      raise SAMPClientError("Client not callable.")

  def register(self):
    """
    Register the client to the SAMP Hub. If the registration fails a L{SAMPClientError}
    is reaised.
    """
    if self.hub.isConnected():

      if self._private_key != None:
        raise SAMPClientError("Client already registered")

      try:
        result = self.hub.register(self.hub.lockfile["samp.secret"])
        if result["samp.self-id"] == "" or result["samp.private-key"] == "":
          raise SAMPClientError("Registation failed. Probably the secret code is wrong.")
        self._public_id = result["samp.self-id"]
        self._private_key = result["samp.private-key"]
        self._hub_id = result["samp.hub-id"]
        if self._callable:
          self._setXmlrpcCallback()
          self._declareSubscriptions()
        if self._metadata != {}:
          self.declareMetadata()
      except SAMPProxyError, err:
        raise SAMPClientError(err.faultString)
      except:
        raise SAMPClientError("Unexpected error: registration failed")

    else:
      raise SAMPClientError("Unable to register to the SAMP Hub. Hub proxy not connected.")

  def unregister(self):
    """
    Unregister the client from the SAMP Hub. If the unregistration fails a L{SAMPClientError}
    is reaised.
    """
    if self.hub.isConnected():

      try:
        self.hub.unregister(self._private_key)
        self._hub_id = None
        self._public_id = None
        self._private_key = None
      except:
        raise SAMPClientError("Unable to unregister from the SAMP Hub.")
    else:
      raise SAMPClientError("Unable to unregister from the SAMP Hub. Hub proxy not connected.")


  def _setXmlrpcCallback(self):
    if self.hub.isConnected() and self._private_key != None:

      try:
        self.hub.setXmlrpcCallback(self._private_key, \
                                   self._xmlrpcAddr)
      except:
        pass

  def _declareSubscriptions(self, subscriptions = None):
    if self.hub.isConnected() and self._private_key != None:

      try:
        mtypes_dict = {}
        # Collect notification mtypes and metadata
        for mtype in self._notification_bindings.keys():
          mtypes_dict[mtype] = copy.deepcopy(self._notification_bindings[mtype][1])

        # Collect notification mtypes and metadata
        for mtype in self._call_bindings.keys():
          mtypes_dict[mtype] = copy.deepcopy(self._call_bindings[mtype][1])

        # Add optional subscription map
        if subscriptions:
          mtypes_dict.update(copy.deepcopy(subscriptions))

        self.hub.declareSubscriptions(self._private_key, mtypes_dict)

      except Exception, ex:
        raise SAMPClientError("Unable to declare subscriptions. Hub unreachable or not connected or client not registered (%s)."%str(ex))
    else:
      raise SAMPClientError("Unable to declare subscriptions. Hub unreachable or not connected or client not registered.")

  def declareMetadata(self, metadata = None):
    """
    Declare the client application metadata supported.

    @param metadata: (optional) dictionary containig the client application metadata
    as defined in the SAMP definition document. If omitted, then none metadata are
    declared.
    @type metadata: dict
    """
    if self.hub.isConnected() and self._private_key != None:

      try:
        if metadata != None:
          self._metadata.update(metadata)

        self.hub.declareMetadata(self._private_key, self._metadata)
      except:
        raise SAMPClientError("Unable to declare metadata. Hub unreachable or not connected or client not registered.")
    else:
      raise SAMPClientError("Unable to declare metadata. Hub unreachable or not connected or client not registered.")

  def getPrivateKey(self):
    """
    Return the client private key used for the Standard Profile communications 
    obtained at registration time (C{samp.private-key}).

    @return: the client private key
    @rtype: string
    """
    return self._private_key

  def getPublicId(self):
    """
    Return public client ID obtained at registration time (C{samp.self-id}).

    @return: the client public ID
    @rtype: string
    """
    return self._public_id

class SAMPHubError(Exception):
  """
    SAMP Hub exceptions.
    """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)


class SAMPClientError(Exception):
  """
    SAMP Client exceptions.
    """
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return repr(self.value)

#: SAMP Proxy Hub exceptions (overwites xmlrpclib.Fault).
SAMPProxyError = xmlrpclib.Fault


class SAMPIntegratedClient(object):
  """
  This class is meant to simplify the client usage providing
  a proxy class that merges the L{SAMPClient} and L{SAMPHubProxy}
  functionalities in a simplified API.
  """
  def __init__(self, name = None, description = None, metadata = None,
               addr = None, port = 0, https = False, key_file=None,
               cert_file=None, cert_reqs=0, ca_certs=None, ssl_version=2,
               callable = True):
    """
    L{SAMPIntegratedClient} constructor.

    @param name: (optional) a string containing the client name
    (corresponding to C{samp.name} metadata keyword).
    @type name: string
    
    @param description: (optional) a string containing the client description
    (corresponding to C{samp.description.text} metadata keyword).
    @type description: string
    
    @param metadata: (optional) a dictionary containing the client 
    application metadata in the standard SAMP format. If present, C{samp.name}
    keyword and C{samp.description.text} keyword are overwritten by the parameters
    C{name} and C{description}.
    @type metadata: dict

    @param addr: (optional) listening address (or IP)
    @type addr: string

    @param port: (optional) the listening XML-RPC server socket port
    @type port: int

    @param https: set the callable client running on a Secure Sockets Layer connection (HTTPS).
    By default SSL is desabled.
    @type https: boolean

    @param key_file: Set the file containing the private key for SSL connections. If the
    certificate file (C{cert_file}) contains the private key, then C{key_file} can be omitted.
    @type key_file: string

    @param cert_file: Specify the file which contains a certificate to be used to identify the
    local side of the secure connection.
    @type cert_file: string

    @param cert_reqs: The parameter C{cert_reqs} specifies whether a certificate is required
    from the Hub side of the connection, and whether it will be validated if provided. It
    must be one of the three values L{ssl.CERT_NONE} (certificates ignored), L{ssl.CERT_OPTIONAL}
    (not required, but validated if provided), or L{ssl.CERT_REQUIRED} (required and validated).
    If the value of this parameter is not L{ssl.CERT_NONE}, then the C{ca_certs} parameter must
    point	to a file of CA certificates.
    @type cert_reqs: int

    @param ca_certs: The C{ca_certs} file contains a set of concatenated \"Certification Authority\" 
    certificates, which are used to validate the certificate passed from the Hub end of the 
    connection.
    @type ca_certs: string

    @param ssl_version: The C{ssl_version} option specifies which version of the SSL protocol to use.
    Typically, the server chooses	a particular protocol version, and the client must adapt to the
    server's choice. Most of the versions are	not interoperable with the other versions. If not
    specified the default SSL version is  L{ssl.PROTOCOL_SSLv23}. This version provides the most
    compatibility with other versions Hub side. Other SSL protocol versions are:                        
    L{ssl.PROTOCOL_SSLv2}, L{ssl.PROTOCOL_SSLv3} and L{ssl.PROTOCOL_TLSv1}.
    @type ssl_version: int

    @param callable: specify whether the client is a callable client or not
    @type callable: boolean

    """

    self.hub = SAMPHubProxy()
    
    self.client = SAMPClient(self.hub, name, description, metadata, addr, port, \
                             https, key_file, cert_file, cert_reqs, \
                             ca_certs, ssl_version, callable)

    
  
    
  def __del__(self):
    try:
      self.disconnect()
    except:
      pass

  # GENERAL
  def isConnected(self):
    """
    Testing method to verify the client connection with a running Hub.

    @return: return True if the client is connected to a Hub, False otherwise
    @rtype: boolean
    """
    return self.hub.isConnected() & self.client.isRunning()

  def connect(self, hub_params = None, user = None, password = None,
              key_file=None, cert_file=None, cert_reqs=0,
              ca_certs=None, ssl_version=1, pool_size=20):
    """
    Connect with the current or specified SAMP Hub, start and register the client.
    If a SAMP Hub is not running or	refuses the connection,	then a L{SAMPHubError} is raised.

    @param hub_params: Optional dictionary containig the lock-file content of the Hub
    with which to connect. This dictionary has the form C{{<token-name>: <token-string>, ...}}.
    @type hub_params: dictionary

    @param user: In case of Basic Authenticated connections, C{user} specifies the user name.
    @type user: string

    @param password: In case of Basic Authenticated connections, C{password} specifies the user
    password.
    @type password: string

    @param key_file: Set the file containing the private key for SSL connections. If the
    certificate file (C{certfile}) contains the private key, then C{keyfile} can be omitted.
    @type key_file: string

    @param cert_file: Specify the file which contains a certificate to be used to identify the
    local side of the secure connection.
    @type cert_file: string

    @param cert_reqs: The parameter C{cert_reqs} specifies whether a certificate is required
    from the server side of the connection, and whether it will be validated if provided. It
    must be one of the three values L{ssl.CERT_NONE} (certificates ignored), L{ssl.CERT_OPTIONAL}
    (not required, but validated if provided), or L{ssl.CERT_REQUIRED} (required and validated).
    If the value of this parameter is not L{ssl.CERT_NONE}, then the C{ca_certs} parameter must
    point	to a file of CA certificates.
    @type cert_reqs: int

    @param ca_certs: The C{ca_certs} file contains a set of concatenated \"Certification Authority\" 
    certificates, which are used to validate the certificate passed from the server end of the 
    connection.
    @type ca_certs: string

    @param ssl_version: The C{ssl_version} option specifies which version of the SSL protocol to use.
    Typically, the server chooses	a particular protocol version, and the client must adapt to the
    server's choice. Most of the versions are	not interoperable with the other versions. If not
    specified the default SSL version is  L{ssl.PROTOCOL_SSLv3}. This version provides the most
    compatibility with other versions server side. Other SSL protocol versions are:                        
    L{ssl.PROTOCOL_SSLv2}, L{ssl.PROTOCOL_SSLv23} and L{ssl.PROTOCOL_TLSv1}.
    @type ssl_version: int

    @param pool_size: The number of socket connections opened to communicate with the Hub.
    @type callable: int
    """
    self.hub.connect(hub_params, user, password, key_file, cert_file, cert_reqs,
                     ca_certs, ssl_version, pool_size)
    self.client.start()
    self.client.register()


  def disconnect(self):
    """
    Unregister the client from the current SAMP Hub, stop the client and disconnect from the Hub.
    """
    cliEx = None
    try:
      self.client.unregister()
    except SAMPClientError, cliEx:
      pass

    if self.client.isRunning():
      self.client.stop()
    self.hub.disconnect()

    if cliEx: raise cliEx


  # HUB
  def ping(self):
    """
    Proxy to C{ping} SAMP Hub method (Standard Profile only)
    """
    return self.hub.ping()

  def declareMetadata(self, metadata):
    """
    Proxy to C{declareMetadata} SAMP Hub method
    """
    return self.client.declareMetadata(metadata)

  def getMetadata(self, client_id):
    """
    Proxy to C{getMetadata} SAMP Hub method
    """
    return self.hub.getMetadata(self.client.getPrivateKey(), client_id)

  def getSubscriptions(self, client_id):
    """
    Proxy to C{getSubscriptions} SAMP Hub method
    """
    return self.hub.getSubscriptions(self.client.getPrivateKey(), client_id)

  def getRegisteredClients(self):
    """
    Proxy to C{getRegisteredClients} SAMP Hub method
    """
    return self.hub.getRegisteredClients(self.client.getPrivateKey())

  def getSubscribedClients(self, mtype):
    """
    Proxy to C{getSubscribedClients} SAMP Hub method
    """
    return self.hub.getSubscribedClients(self.client.getPrivateKey(), mtype)

  def _format_easy_msg(self, mtype, params):

    msg = {}

    if params.has_key("extra_kws"):
      extra = params["extra_kws"]
      del(params["extra_kws"])
      msg = {"samp.mtype": mtype, "samp.params": params}
      msg.update(extra)
    else:
      msg = {"samp.mtype": mtype, "samp.params": params}

    return msg


  def notify(self, recipient_id, message):
    """
    Proxy to C{notify} SAMP Hub method
    """
    return self.hub.notify(self.client.getPrivateKey(), recipient_id, message)

  def enotify(self, recipient_id, mtype, **params):
    """
    Easy C{notify}. It is a proxy to L{notify} method that allows to
    send the notification message in a simplified way. Example:

    >>> import sampy
    >>> cli = sampy.SAMPIntegratedClient()
    >>> ...
    >>> cli.enotify("samp.msg.progress", msgid = "xyz", txt = "initialization", \\
    >>>             percent = "10", extra_kws = {"my.extra.info": "just an example"})

    Note that reserved C{extra_kws} keyword is a dictionary with the special meaning of 
    being used to add extra keywords, in addition to the standard C{samp.mtype}
    and C{samp.params}, to the message sent.

    @param recipient_id: the recipient ID
    @type recipient_id: string

    @param mtype: the MType to be notified
    @type mtype: string

    @param params: variable keyword set which contains the list of parameters for 
                   the specified MType
    @type params: dictionary or set of keywords

    """
    return self.notify(recipient_id, self._format_easy_msg(mtype, params))

  def notifyAll(self, message):
    """
    Proxy to C{notifyAll} SAMP Hub method
    """
    return self.hub.notifyAll(self.client.getPrivateKey(), message)

  def enotifyAll(self, mtype, **params):
    """
    Easy C{notify}. It is a proxy to L{notifyAll} method that allows to
    send the notification message in a simplified way. Example:

    >>> import sampy
    >>> cli = sampy.SAMPIntegratedClient()
    >>> ...
    >>> cli.enotifyAll("samp.msg.progress", txt = "initialization", \\
    >>>                percent = "10", extra_kws = {"my.extra.info": "just an example"})

    Note that reserved C{extra_kws} keyword is a dictionary with the special meaning of 
    being used to add extra keywords, in addition to the standard C{samp.mtype}
    and C{samp.params}, to the message sent.

    @param mtype: the MType to be notified
    @type mtype: string

    @param params: variable keyword set which contains the list of parameters for 
                   the specified MType
    @type params: dictionary or set of keywords

    """
    return self.notifyAll(self._format_easy_msg(mtype, params))


  def call(self, recipient_id, msg_tag, message):
    """
    Proxy to C{call} SAMP Hub method
    """
    return self.hub.call(self.client.getPrivateKey(), recipient_id, msg_tag, message)

  def ecall(self, recipient_id, msg_tag, mtype, **params):
    """
    Easy C{call}. It is a proxy to L{call} method that allows to
    send a call message in a simplified way. Example:

    >>> import sampy
    >>> cli = sampy.SAMPIntegratedClient()
    >>> ...
    >>> msgid = cli.ecall("abc", "xyz", "samp.msg.progress", txt = "initialization", \\
    >>>                   percent = "10", extra_kws = {"my.extra.info": "just an example"})

    Note that reserved C{extra_kws} keyword is a dictionary with the special meaning of 
    being used to add extra keywords, in addition to the standard C{samp.mtype}
    and C{samp.params}, to the message sent.

    @param recipient_id: the recipient ID
    @type recipient_id: string

    @param msg_tag: the message tag to use
    @type msg_tag: string

    @param mtype: the MType to be sent
    @type mtype: string

    @param params: variable keyword set which contains the list of parameters for 
                   the specified MType
    @type params: dictionary or set of keywords
    """

    return self.call(recipient_id, msg_tag, self._format_easy_msg(mtype, params))


  def callAll(self, msg_tag, message):
    """
    Proxy to C{callAll} SAMP Hub method
    """
    return self.hub.callAll(self.client.getPrivateKey(), msg_tag, message)

  def ecallAll(self, msg_tag, mtype, **params):
    """
    Easy C{callAll}. It is a proxy to L{callAll} method that allows to
    send the call message in a simplified way. Example:

    >>> import sampy
    >>> cli = sampy.SAMPIntegratedClient()
    >>> ...
    >>> msgid = cli.ecallAll("xyz", "samp.msg.progress", txt = "initialization", \\
    >>>                      percent = "10", extra_kws = {"my.extra.info": "just an example"})

    Note that reserved C{extra_kws} keyword is a dictionary with the special meaning of 
    being used to add extra keywords, in addition to the standard C{samp.mtype}
    and C{samp.params}, to the message sent.

    @param msg_tag: the message tag to use
    @type msg_tag: string

    @param mtype: the MType to be sent
    @type mtype: string

    @param params: variable keyword set which contains the list of parameters for 
                   the specified MType
    @type params: dictionary or set of keywords
    """
    self.callAll(msg_tag, self._format_easy_msg(mtype, params))


  def callAndWait(self, recipient_id, message, timeout):
    """
    Proxy to C{callAndWait} SAMP Hub method. If timeout expires a 
    L{SAMPProxyError} instance is raised.
    """
    return self.hub.callAndWait(self.client.getPrivateKey(), recipient_id, message, timeout)

  def ecallAndWait(self, recipient_id, mtype, timeout, **params):
    """
    Easy C{callAndWait}. It is a proxy to L{callAll} method that allows to
    send the call message in a simplified way. Example:

    >>> import sampy
    >>> cli = sampy.SAMPIntegratedClient()
    >>> ...
    >>> cli.ecallAndWait("xyz", "samp.msg.progress", "5", txt = "initialization", \\
    >>>                  percent = "10", extra_kws = {"my.extra.info": "just an example"})

    Note that reserved C{extra_kws} keyword is a dictionary with the special meaning of 
    being used to add extra keywords, in addition to the standard C{samp.mtype}
    and C{samp.params}, to the message sent.

    @param recipient_id: the recipient ID
    @type recipient_id: string

    @param mtype: the MType to be sent
    @type mtype: string

    @param timeout: the call timeout in seconds
    @type timeout: string

    @param params: variable keyword set which contains the list of parameters for 
                   the specified MType
    @type params: dictionary or set of keywords
    """
    return self.callAndWait(recipient_id, self._format_easy_msg(mtype, params), timeout)


  def reply(self, msg_id, response):
    """
    Proxy to C{reply} SAMP Hub method
    """
    return self.hub.reply(self.client.getPrivateKey(), msg_id, response)

  def _format_easy_response(self, status, result, error):

    msg = {"samp.status": status}
    if result != None:
      msg.update({"samp.result": result})
    if error != None:
      msg.update({"samp.error": error})

    return msg

  def ereply(self, msg_id, status, result = None, error = None):
    """
    Easy C{reply}. It is a proxy to L{callAll} method that allows to
    send a reply message in a simplified way. Example:

    >>> import sampy
    >>> cli = sampy.SAMPIntegratedClient()
    >>> ...
    >>> cli.ereply("abd", sampy.SAMP_STATUS_ERROR, result = {}, error = {"samp.errortxt": "Test error message"})

    @param msg_id: the message ID to which reply
    @type msg_id: string

    @param status: the content of the C{samp.status} response keyword
    @type status: string

    @param result: the content of the C{samp.result} response keyword
    @type result: dictionary

    @param error: the content of the C{samp.error} response keyword
    @type error: dictionary
    """
    return self.reply(msg_id, self._format_easy_response(status, result, error))




  # CLIENT
  def receiveNotification(self, private_key, sender_id, message):
    """
    Standard callable client C{receiveNotification} method. This method is
    automatically handled when L{bindReceiveNotification} method is used to bind
    distinct operations to MTypes. In case of a customized callable client
    implementation that inherits from L{SAMPClient} class this method should
    be overwritten. ATTENTION: when overwritten, this method must always return
    a string result (even empty).

    @param private_key: the client private key.
    @type private_key: str

    @param sender_id: the sender public ID.
    @type sender_id: str

    @param message: the message received.
    @type message: dict

    @return: any confirmation string.
    @rtype: str
    """
    return self.client.receiveNotification(private_key, sender_id, message)

  def receiveCall(self, private_key, sender_id, msg_id, message):
    """
    Standard callable client C{receiveCall} method. This method is
    automatically handled when L{bindReceiveCall} method is used to bind
    distinct operations to MTypes. In case of a customized callable client
    implementation that inherits from L{SAMPClient} class this method should
    be overwritten. ATTENTION: when overwritten, this method must always return
    a string result (even empty).

    @param private_key: the client private key.
    @type private_key: str

    @param sender_id: the sender public ID.
    @type sender_id: str

    @param msg_id: the message ID received.
    @type msg_id: str

    @param message: the message received.
    @type message: dict

    @return: any confirmation string.
    @rtype: str
    """
    return self.client.receiveCall(private_key, sender_id, msg_id, message)

  def receiveResponse(self, private_key, responder_id, msg_tag, response):
    """
    Standard callable client C{receiveResponse} method. This method is
    automatically handled when L{bindReceiveResponse} method is used to bind
    distinct operations to MTypes. In case of a customized callable client
    implementation that inherits from L{SAMPClient} class this method should
    be overwritten. ATTENTION: when overwritten, this method must always return
    a string result (even empty).

    @param private_key: the client private key.
    @type private_key: str

    @param responder_id: the responder public ID.
    @type responder_id: str

    @param msg_tag: the response message tag.
    @type msg_tag: str

    @param response: the response received.
    @type response: dict

    @return: any confirmation string.
    @rtype: str
    """
    return self.client.receiveResponse(private_key, responder_id, msg_tag, response)



  def bindReceiveMessage(self, mtype, function, declare = True, metadata = None):
    """Bind a specific MType to a function or class method, being intended for
    a call or a notification.

    The function must be of the form:
    C{def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)}

    where C{private_key} is the client private-key, C{sender_id} argument is the
    notification sender ID, C{msg_id} is the Hub message-id (calls only, otherwise is None),
    C{mtype} is the message MType, C{params} is the message parameter set (content of
    "samp.params") and C{extra} is a dictionary containing any extra message map entry.
    The client is automatically declared subscribed to the MType by default.

    @param mtype: the MType to be catched.
    @type mtype: str

    @param function: the application function to be used when C{mtype} is received.
    @type function: function or class method

    @param declare: specify whether the client must be automatically declared as
    subscribed to the MType (see also L{declareSubscriptions}).
    @type declare: boolean

    @param metadata: an optional map containing additional metadata to declare associated
    with the MType subscribed to (see also L{declareSubscriptions}).
    @type metadata: dict
    """
    self.client.bindReceiveMessage(mtype, function, declare = True, metadata = None)

  def bindReceiveNotification(self, mtype, function, declare = True, metadata = None):
    """
    Bind a specific MType notification to a function or class method.
    The function must be of the form:

    C{def my_function_or_method(<self,> private_key, sender_id, mtype, params, extra)}

    where C{private_key} is the client private-key, C{sender_id} argument is 
    the notification sender ID, C{mtype} is the message MType, C{params} is 
    the notified message parameter set (content of "samp.params") and C{extra} is a 
    dictionary containing any extra message map entry. The client is
    automatically declared subscribed to the MType by default.

    @param mtype: the MType to be catched.
    @type mtype: str

    @param function: the application function to be used when C{mtype} is received.
    @type function: function or class method

    @param declare: specify whether the client must be automatically declared as
    subscribed to the MType (see alse L{declareSubscriptions}).
    @type declare: boolean

    @param metadata: an optional map containing additional metadata to declare associated
    with the MType subscribed to (see also L{declareSubscriptions}).
    @type metadata: dict
    """

    self.client.bindReceiveNotification(mtype, function, declare, metadata)

  def bindReceiveCall(self, mtype, function, declare = True, metadata = None):
    """
    Bind a specific MType call to a function or class method.
    The function must be of the form:

    C{def my_function_or_method(<self,> private_key, sender_id, msg_id, mtype, params, extra)}

    where C{private_key} is the client private-key, C{sender_id} argument is the
    notification sender ID, C{msg_id} is the Hub message-id, C{mtype} is the message MType, 
    C{params} is the message parameter set (content of "samp.params") and C{extra} is a 
    dictionary containing any extra message map entry. The client is
    automatically declared subscribed to the MType by default.

    @param mtype: the MType to be catched.
    @type mtype: str

    @param function: the application function to be used when C{mtype} is received.
    @type function: function or class method

    @param declare: specify whether the client must be automatically declared as
    subscribed to the MType (see also L{declareSubscriptions}).
    @type declare: boolean

    @param metadata: an optional map containing additional metadata to declare associated
    with the MType subscribed to (see also L{declareSubscriptions}).
    @type metadata: dict
    """
    self.client.bindReceiveCall(mtype, function, declare, metadata)

  def bindReceiveResponse(self, msg_tag, function):
    """
    Bind a specific msg-tag response to a function or class method.
    The function must be of the form:

    C{def my_function_or_method(<self,> private_key, responder_id, msg_tag, response)}

    where C{private_key} is the client private-key, C{responder_id} argument is the message
    responder ID, C{msg_tag} is the message-tag provided at call time and C{response} is the
    response received.

    @param msg_tag: the message-tag to be catched.
    @type msg_tag: str

    @param function: the application function to be used when C{msg_tag} is received.
    @type function: function or class method
    """
    self.client.bindReceiveResponse(msg_tag, function)

  def unbindReceiveNotification(self, mtype, declare = True):
    """
    Remove from the notifications binding table the specified MType and unsubscribe
    the client from it (if required).

    @param mtype: the MType to be removed
    @type mtype: str

    @param declare: specify whether the client must be automatically declared as
    unsubscribed from the MType (see alse L{declareSubscriptions}).
    @type declare: boolean

    """
    self.client.unbindReceiveNotification(mtype, declare)

  def unbindReceiveCall(self, mtype, declare = True):
    """
    Remove from the calls binding table the specified MType and unsubscribe
    the client from it (if required).

    @param mtype: the MType to be removed
    @type mtype: str

    @param declare: specify whether the client must be automatically declared as
    unsubscribed from the MType (see alse L{declareSubscriptions}).
    @type declare: boolean
    """
    self.client.unbindReceiveCall(mtype, declare)

  def unbindReceiveResponse(self, msg_tag):
    """
    Remove from the responses binding table the specified message-tag.

    @param msg_tag: the message-tag to be removed
    @type msg_tag: str
    """
    self.client.unbindReceiveResponse(msg_tag)


  def declareSubscriptions(self, subscriptions = None):
    """
    Declares the MTypes the client wishes to subscribe to, implicitly defined
    with the MType binding methods L{bindReceiveNotification} and L{bindReceiveCall}.
    An optional C{subscriptions} map can be added to the final map passed to 
    the L{SAMPHubProxy.declareSubscriptions} operation.

    @param subscriptions: an optional map containing the list of MTypes to subscribe to,
    with the same format of the C{subscriptions} map passed to the
    L{SAMPHubProxy.declareSubscriptions} operation.
    @type subscriptions: dict
    """
    self.client.declareSubscriptions(subscriptions)

  def getPrivateKey(self):
    """
    Return the client private key used for the Standard Profile communications 
    obtained at registration time (C{samp.private-key}).

    @return: the client private key
    @rtype: string
    """
    return self.client.getPrivateKey()

  def getPublicId(self):
    """
    Return public client ID obtained at registration time (C{samp.self-id}).

    @return: the client public ID
    @rtype: string
    """
    return self.client.getPublicId()



def main(timeout=0):

  import signal
  from optparse import OptionParser, OptionGroup
  parser = OptionParser(version="%prog " + __release__)

  parser.disable_interspersed_args()

  parser.add_option("-k", "--secret", dest="secret", metavar="CODE",
                    help="custom secret code.")

  parser.add_option("-d", "--addr", dest="addr", metavar="ADDR",
                    help="listening address (or IP).")

  parser.add_option("-p", "--port", dest="port", metavar="PORT", type="int",
                    help="listening port number.")

  parser.add_option("-f", "--lockfile", dest="lockfile", metavar="FILE",
                    help="custom lockfile.")

  parser.add_option("-w", "--no-web-profile", dest="web_profile", action = "store_false",
                    help = "run the Hub disabling the Web Profile.")

  parser.add_option("-P", "--pool-size", dest="pool_size", metavar="SIZE", type="int",
                    help = "the socket connections pool size.")
  
  timeout_group = OptionGroup(parser, "Timeout group",
                              "Special options to setup hub and client timeouts." \
                              "It contains a set of special options that allows to set up the Hub and " \
                              "clients inactivity timeouts, that is the Hub or client inactivity time " \
                              "interval after which the Hub shuts down or unregisters the client. " \
                              "Notification of samp.hub.disconnect MType is sent to the clients " \
                              "forcibly unregistered for timeout expiration.")

  timeout_group.add_option("-t", "--timeout", dest="timeout", metavar="SECONDS",
                           help="set the Hub inactivity timeout in SECONDS. By default it " \
                           "is set to 0, that is the Hub never expires.", type= "int")


  timeout_group.add_option("-c", "--client-timeout", dest="client_timeout", metavar="SECONDS",
                           help="set the client inactivity timeout in SECONDS. By default it "\
                           "is set to 0, that is the client never expires.", type= "int")

  parser.add_option_group(timeout_group)

  log_group = OptionGroup(parser, "Logging options",
                          "Additional options which allow to customize the logging output. By " \
                          "default SAMPy Hub uses the standard output and standard error " \
                          "devices to print out INFO level logging messages. Using the options " \
                          "here below it is possible to modify the logging level and also " \
                          "specify the output files where redirect the logging messages.")

  log_group.add_option("-L", "--log-level", dest="loglevel", metavar="LEVEL",
                       help="set the Hub instance log level (OFF, ERROR, WARNING, INFO, DEBUG).",
                       type = "choice", choices = ["OFF", "ERROR", "WARNING", "INFO", "DEBUG"])


  log_group.add_option("-O", "--log-output", dest="logout", metavar="FILE",
                       help="set the output file for INFO and DEBUG messages.")


  log_group.add_option("-E", "--log-error", dest="logerr", metavar="FILE",
                       help="set the output file for WARNING and ERROR messages.")


  parser.add_option_group(log_group)

  adv_group = OptionGroup(parser, "Advanced group",
                          "Advanced options addressed to facilitate administrative tasks and " \
                          "allow new non-standard Hub behaviors. In particular the --label " \
                          "options is used to assign a value to hub.label token and is used to " \
                          "assign a name to the Hub instance. Option --auth-file sets the " \
                          "authentication file enabling the access control through Basic " \
                          "Authentication. The authentication file is a Berkeley DB file " \
                          "in Hash format containing a set of <user name>=md5(<password>)<group 1>,<group 2>,<group 3>,... " \
                          "key=value pairs. Options --owner and --group are additional tokens " \
                          "(hub.owner.name and hub.owner.group respectively) available for " \
                          "possible administrative tasks and usable in conjunction with --restrict " \
                          "option for authentication tasks. --restrict option allows to restrict "\
                          "the Hub access only to OWNER or to a certain owner GROUP. Access restriction " \
                          "is actually enabled only if an authentication file is provided. " \
                          "--admin option defines the name of the administrator user in case " \
                          "of restricted access. The administrator user can always access the " \
                          "hub instance even if it is running with the OWNER restricion policy. " \
                          "The administrator must be declared in the authentication file. " \
                          "The very special --multi option allows to start a Hub in multi-instance mode. " \
                          "Multi-instance mode is a non-standard Hub behavior that enables " \
                          "multiple contemporaneous running Hubs. Multi-instance hubs place " \
                          "their non-standard lock-files within the <home directory>/.samp-1 " \
                          "directory naming them making use of the format: " \
                          "samp-hub-<PID>-<ID>, where PID is the Hub process ID while ID is an " \
                          "internal ID (integer).")

  adv_group.add_option("-l", "--label", dest = "label", metavar = "LABEL",
                       help = "assign a LABEL to the Hub.")	

  adv_group.add_option("-o", "--owner", dest="owner", metavar="OWNER",
                       help="assign a OWNER to the Hub.")


  adv_group.add_option("-g", "--group", dest = "owner_group", metavar = "GROUP",
                       help = "assign a GROUP to the Hub.")

  if BDB_SUPPORT:
    adv_group.add_option("-a", "--auth-file", dest = "auth_file", metavar = "FILE",
                         help = "filename for Basic Authentication.")	

    adv_group.add_option("-r", "--restrict", dest = "access_restrict", metavar = "LEVEL",
                         help = "set the access restriction to OWNER or GROUP.",
                         type = "choice", choices = ["OWNER", "GROUP"])

    adv_group.add_option("-A", "--admin", dest = "admin", metavar = "NAME",
                         help = "set the name of the administrator user.")

  adv_group.add_option("-m", "--multi", dest = "mode",
                       help = "run the Hub in multi-instance mode generating a custom " \
                       "lockfile with a random name.",
                       action = "store_const", const = SAMP_HUB_MULTIPLE_INSTANCE)


  parser.add_option_group(adv_group)


  if SSL_SUPPORT:

    ssl_group = OptionGroup(parser, "SSL group", "Additional options to launch " \
                            "the Hub instance using the Secure Sockets Layer (HTTPS). " \
                            "The --key-file and --cert-file parameters specify optional " \
                            "files which contain a certificate to be used to identify " \
                            "the local side of the connection. " \
                            "Often the private key is stored in the same file as the " \
                            "certificate; in this case, only the --cert-file parameter need " \
                            "be passed. If the private key is stored in a separate file, " \
                            "both parameters must be used. If the private key is stored " \
                            "in the certificate file, it should come before the first certificate " \
                            "in the certificate chain.")

    ssl_group.add_option("-s", "--https", dest = "https", action = "store_true",
                         help = "run the Hub using the Secure Sockets Layer.")

    ssl_group.add_option("-C", "--cert-file", dest = "certfile", metavar = "FILE",
                         help = "set the certificate file.")

    ssl_group.add_option("-K", "--key-file", dest = "keyfile", metavar = "FILE",
                         help = "set the key file. By default this option is ignored, " \
                         "assuming that the private key is stored in the certificate file.")

    ssl_group.add_option("--cert-reqs", dest = "cert_reqs", metavar = "STRING",
                         help = "this option specifies whether a certificate " \
                         "is required from the client side of the connection, and whether " \
                         "it will be validated if provided. It must be one of the three " \
                         "values NONE (certificates ignored, default), OPTIONAL (not " \
                         "required, but validated if provided), or REQUIRED (required " \
                         "and validated). If the value of this option is not NONE, " \
                         "then the --ca-certs option must point to a file of CA certificates.",
                         type = "choice", choices = ["NONE", "OPTIONAL", "REQUIRED"])

    ssl_group.add_option("--ca-certs", dest = "ca_certs", metavar = "FILE",
                         help = "the --ca-certs file contains a set of concatenated " \
                         "\"certification authority\" certificates, which are used to " \
                         "validate certificates passed from the client end of the " \
                         "connection.")

    ssl_group.add_option("--ssl-version", dest = "ssl_version", metavar = "STRING",
                         help = "the --ssl-version option specifies which version of the " \
                         "SSL protocol to use. Typically, the server chooses a particular " \
                         "protocol version, and the client must adapt to the server's choice. " \
                         "Most of the versions are not interoperable with the other versions. " \
                         "If not specified the default SSL version is SSLv23. This version " \
                         "provides the most compatibility with other versions client side. " \
                         "Other SSL protocol versions are: SSLv2, SSLv3 and TLSv1.",
                         type = "choice", choices = ["SSLv23", "SSLv2", "SSLv3", "TLSv1"])

    parser.add_option_group(ssl_group)

    parser.set_defaults(https = False)
    parser.set_defaults(certfile = None)
    parser.set_defaults(keyfile = None)
    parser.set_defaults(cert_reqs = "NONE")
    parser.set_defaults(ca_certs = None)
    parser.set_defaults(ssl_version = "SSLv23")

  parser.set_defaults(web_profile = True)
  parser.set_defaults(pool_size = 20)
  parser.set_defaults(timeout = 0)
  parser.set_defaults(client_timeout = 0)
  parser.set_defaults(loglevel = "INFO")
  parser.set_defaults(logout = "")
  parser.set_defaults(logerr = "")
  parser.set_defaults(label = "")
  parser.set_defaults(owner = "")
  parser.set_defaults(owner_group = "")
  parser.set_defaults(admin = "admin")

  if BDB_SUPPORT:
    parser.set_defaults(auth_file = None)
    parser.set_defaults(access_restrict = None)

  parser.set_defaults(mode = SAMP_HUB_SINGLE_INSTANCE)		

  (options, args) = parser.parse_args()

  try:


    if SSL_SUPPORT:

      # Set ssl options properly

      if options.cert_reqs == "OPTIONAL":
        options.cert_reqs = ssl.CERT_OPTIONAL
      elif options.cert_reqs == "REQUIRED":
        options.cert_reqs = ssl.CERT_REQUIRED
      else:
        options.cert_reqs = ssl.CERT_NONE

      if options.ssl_version == "SSLv2":
        options.ssl_version = ssl.PROTOCOL_SSLv2
      elif options.ssl_version == "SSLv3":
        options.ssl_version = ssl.PROTOCOL_SSLv3
      elif options.ssl_version == "TLSv1":
        options.ssl_version = ssl.PROTOCOL_TLSv1
      else:
        options.ssl_version = ssl.PROTOCOL_SSLv23

    level = SAMPLog.INFO
    stdout = sys.stdout
    stderr = sys.stderr

    if options.loglevel in ("OFF", "ERROR", "WARNING", "DEBUG"):
      if options.loglevel == "OFF":
        level = SAMPLog.OFF
      if options.loglevel == "ERROR":
        level = SAMPLog.ERROR
      if options.loglevel == "WARNING":
        level = SAMPLog.WARNING
      if options.loglevel == "DEBUG":
        level = SAMPLog.DEBUG

    if options.logout != "":
      stdout = open(options.logout, "a")
    if options.logerr != "":
      stderr = open(options.logerr, "a")

    log = SAMPLog(level = level, stdout = stdout, stderr = stderr)

    if BDB_SUPPORT:

      if options.access_restrict != None:
        if options.access_restrict == SAMP_RESTRICT_OWNER and options.owner == "":
          log.warning("The access cannot be restricted to the owner if the owner "
                      "name is not specified!")
          sys.exit()

        if options.access_restrict == SAMP_RESTRICT_GROUP and options.owner_group == "":
          log.warning("The access cannot be restricted to the owner group if the owner "
                      "group name is not specified!")
          sys.exit()

    args = copy.deepcopy(options.__dict__)
    del(args["loglevel"])
    del(args["logout"])
    del(args["logerr"])
    args["log"] = log

    hub = SAMPHubServer(**args)
    hub.start(False)

    if not timeout:
      while hub.isRunning(): time.sleep(0.01)
    else:
      time.sleep(timeout)
      hub.stop()

  except KeyboardInterrupt:
    hub.stop()
  except IOError, (errno, strerror):
    print "[SAMP] Error: I/O error(%s): %s" % (errno, strerror)
  except SAMPHubError:
    pass
  except SystemExit:
    pass
  except:
    print "[SAMP] Error: Unexpected error:", sys.exc_info()


if __name__ == "__main__":

  main()
