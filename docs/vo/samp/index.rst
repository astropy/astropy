.. _vo-samp:

***************************************************************
SAMP (Simple Application Messaging Protocol (`astropy.vo.samp`)
***************************************************************

Introduction
============

`astropy.vo.samp` is an IVOA SAMP (Simple Application Messaging Protocol) messaging system implementation in Python.

Messaging Protocol version 1.3.

It provides classes to easily:

1. instantiate one or multiple Hubs;
2. interface an application or script to a running Hub;
3. create and manage a SAMP client.

`astropy.vo.samp` provides also a stand-alone program ``sampy`` capable to
instantiate a persistent Hub. In order to have a full description of
``sampy`` stand-alone program options, type the command::

   $ sampy -h


.. _IVOA Simple Application Messaging Protocol: http://www.ivoa.net/Documents/latest/SAMP.html


Getting Started
===============


Using `astropy.vo.samp`
=======================


Examples
^^^^^^^^

This module contains classes to create a SAMP Hub and/or interface 
an application to a running SAMP Hub (Standard Profile).

SAMPy is very simple to use for the creation of a running SAMP Hub. The only thing
you have to do is to create a ``SAMPHubServer`` instance and start it:

>>> from astropy.vo.samp import *
>>> hub = SAMPHubServer()
>>> hub.start()

``sampy.py`` module file can also be used as a stand-alone program to launch
new SAMP Hub instances. Just run it::

   $ ./sampy.py

To have a more exhaustive documentation on ``sampy.py`` as stand-alone program,
please read the ``README`` file distributed with SAMPy.

SAMPy can also be used to easily interface your Python script (or application) to a
running SAMP Hub and create a callable Client thanks to the ``SAMPIntegratedClient``
class which integrates the funcionalities of the more specific classes ``SAMPClient``
and ``SAMPHubProxy``.

>>> from astropy.vo.samp import *
>>> import time
>>>>>>>>>>>>>>>
>>> # Create a client
... cli1 = SAMPIntegratedClient(name = "Client 1", description = "Test Client 1",
...                             metadata = {"cli1.version":"0.01"})
>>> # Create another client (name and description included in metadata argument)
... cli2 = SAMPIntegratedClient(metadata = {"samp.name":"Client 2",
...                                         "samp.description.text":"Test Client 2",
...                                         "cli2.version":"0.25"})
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>> # Connect them
... cli1.connect()
>>> cli2.connect()
>>>>>>>>>>>>>>>>>>
>>>
>>> print("CLI1", cli1.getPrivateKey(), cli1.getPublicId())
>>> print("CLI2", cli2.getPrivateKey(), cli2.getPublicId())
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>
>>> # Function called when a notification is received
... def test_receive_notification(private_key, sender_id, mtype, params, extra):
...   print("Notification:", private_key, sender_id, mtype, params, extra)
..........................................................................
>>> # Function called when a call is received
... def test_receive_call(private_key, sender_id, msg_id, mtype, params, extra):
...   print("Call:", private_key, sender_id, msg_id, mtype, params, extra)
...   cli1.ereply(msg_id, SAMP_STATUS_OK, result = {"txt": "printed"})
......................................................................
>>> # Function called when a response is received
... def test_receive_response(private_key, sender_id, msg_id, response):
...   print("Response:", private_key, sender_id, msg_id, response)
..................................................................
>>> # Subscribe Client 1 to "samp.*" and "samp.app.*" MType and bind it to
... # the related functions
... cli1.bindReceiveNotification("samp.app.*", test_receive_notification)
>>> cli1.bindReceiveCall("samp.app.*", test_receive_call)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>> # Bind Client 2 message-tags received to suitable functions
... cli2.bindReceiveResponse("my-dummy-print", test_receive_response)
>>> cli2.bindReceiveResponse("my-dummy-print-specific", test_receive_response)
>>> # Client 2 notifies to All "samp.app.echo" MType using myhub
... cli2.enotifyAll("samp.app.echo", txt = "Hello world!")
['cli#2']
>>> Notification: 0d7f4500225981c104a197c7666a8e4e cli#2 samp.app.echo {'txt': 'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>> print(cli2.getSubscribedClients("samp.app.echo"))
{'cli#2': {}}
>>> # Client 2 calls to All "samp.app.echo" MType using "my-dummy-print"
... # as message-tag
... print(cli2.callAll("my-dummy-print",
...                    {"samp.mtype": "samp.app.echo",
...                     "samp.params": {"txt": "Hello world!"}}))
{'cli#1': 'msg#1;;cli#hub;;cli#2;;my-dummy-print'}
>>> Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2 msg#1;;cli#hub;;cli#2;;my-dummy-print samp.app.echo {'txt': 'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>> # Client 2 calls "samp.app.echo" MType on Client 1 tagging it as
... # "my-dummy-print-specific"
... try:
...   print(cli2.call(cli1.getPublicId(),
...                   "my-dummy-print-specific",
...                   {"samp.mtype": "samp.app.echo",
...                    "samp.params": {"txt": "Hello Cli 1!"}}))
... except SAMPProxyError, e:
...   print("Error (%d): %s" % (e.faultCode, e.faultString))
............................................................
msg#2;;cli#hub;;cli#2;;my-dummy-print-specific
>>> Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2 msg#2;;cli#hub;;cli#2;;my-dummy-print-specific samp.app.echo {'txt': 'Hello Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print-specific {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>
>>>
>>>
>>> # Function called to test synchronous calls
... def test_receive_sync_call(private_key, sender_id, msg_id, mtype, params, extra):
...   print("SYNC Call:", sender_id, msg_id, mtype, params, extra)
...   time.sleep(2)
...   cli1.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
...                       "samp.result": {"txt": "printed sync"}})
..................................................................
>>>
>>> # Bind test MType for sync calls
... cli1.bindReceiveCall("samp.test", test_receive_sync_call)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>>
>>>
>>>
>>> try:
...   # Sync call
...   print(cli2.callAndWait(cli1.getPublicId(),
...                          {"samp.mtype": "samp.test",
...                           "samp.params": {"txt": "Hello SYNCRO Cli 1!"}},
...                          "10"))
... except SAMPProxyError, e:
...   # If timeout expires than a SAMPProxyError is returned
...   print("Error (%s): %s" % (e.faultCode, e.faultString))
............................................................
SYNC Call: cli#2 msg#3;;cli#hub;;cli#2;;sampy::sync::call samp.test {'txt': 'Hello SYNCRO Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
{'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed sync'}}
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
>>> cli1.disconnect()
>>> cli2.disconnect()


Reference/API
=============

.. automodapi:: astropy.vo.samp


Acknowledgments and Licenses
============================

