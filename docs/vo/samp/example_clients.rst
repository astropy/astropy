.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example_clients:

*****************************************************
Using ``SAMPIntegratedClient`` objects to communicate
*****************************************************

``astropy.vo.samp`` can also be used to easily interface your Python script (or application) to a
running SAMP Hub and create a callable client thanks to the |SAMPIntegratedClient|
class which integrates the funcionalities of the more specific classes |SAMPClient|
and |SAMPHubProxy|.

First start a SAMP hub::

   >>> from astropy.vo import samp
   >>> hub = samp.SAMPHubServer()
   >>> hub.start()

Next create two clients and connect them::

   >>> client1 = samp.SAMPIntegratedClient(name = "Client 1", description = "Test Client 1",
   ...                                     metadata = {"client1.version":"0.01"})
   >>> client2 = samp.SAMPIntegratedClient(metadata = {"samp.name":"Client 2",
   ...                                                 "samp.description.text":"Test Client 2",
   ...                                                 "client2.version":"0.25"})
   >>> client1.connect()
   >>> client2.connect()
   >>> print "client1", client1.getPrivateKey(), client1.getPublicId()
   >>> print "client2", client2.getPrivateKey(), client2.getPublicId()

Define functions to call when receiving a notification, call or response::

   >>> def test_receive_notification(private_key, sender_id, mtype, params, extra):
   ...     print("Notification:", private_key, sender_id, mtype, params, extra)
   >>> def test_receive_call(private_key, sender_id, msg_id, mtype, params, extra):
   ...     print("Call:", private_key, sender_id, msg_id, mtype, params, extra)
   ...     client1.ereply(msg_id, SAMP_STATUS_OK, result = {"txt": "printed"})
   >>> def test_receive_response(private_key, sender_id, msg_id, response):
   ...     print("Response:", private_key, sender_id, msg_id, response)

Subscribe client 1 to ``"samp.*"`` and ``"samp.app.*"`` MTypes and bind them to the related functions::

   >>> client1.bindReceiveNotification("samp.app.*", test_receive_notification)
   >>> client1.bindReceiveCall("samp.app.*", test_receive_call)


Bind client 2 message-tags received to suitable functions::

   >>> client2.bindReceiveResponse("my-dummy-print", test_receive_response)
   >>> client2.bindReceiveResponse("my-dummy-print-specific", test_receive_response)


Client 2 notifies to All "samp.app.echo" MType using the hub::

   >>> client2.enotifyAll("samp.app.echo", txt="Hello world!")
   ['cli#2']
   Notification: 0d7f4500225981c104a197c7666a8e4e cli#2 samp.app.echo {'txt': 'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
   >>> print(client2.getSubscribedClients("samp.app.echo"))
   {'cli#2': {}}

Client 2 calls to all ``"samp.app.echo"`` MType using ``"my-dummy-print"`` as message-tag::

   >>> print(client2.callAll("my-dummy-print",
   ...                       {"samp.mtype": "samp.app.echo",
   ...                        "samp.params": {"txt": "Hello world!"}}))
   {'cli#1': 'msg#1;;cli#hub;;cli#2;;my-dummy-print'}
   Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2 msg#1;;cli#hub;;cli#2;;my-dummy-print samp.app.echo {'txt': 'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
   Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}


Client 2 calls ``"samp.app.echo"`` MType on client 1 tagging it as ``"my-dummy-print-specific"``::

   >>> try:
   ...     print(client2.call(client1.getPublicId(),
   ...                        "my-dummy-print-specific",
   ...                        {"samp.mtype": "samp.app.echo",
   ...                         "samp.params": {"txt": "Hello client 1!"}}))
   ... except SAMPProxyError as e:
   ...     print("Error ({0}): {1}".format(e.faultCode, e.faultString))
   msg#2;;cli#hub;;cli#2;;my-dummy-print-specific
   Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2 msg#2;;cli#hub;;cli#2;;my-dummy-print-specific samp.app.echo {'txt': 'Hello Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
   Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print-specific {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}


Function called to test synchronous calls::

   >>> def test_receive_sync_call(private_key, sender_id, msg_id, mtype, params, extra):
   ...     import time
   ...     print("SYNC Call:", sender_id, msg_id, mtype, params, extra)
   ...     time.sleep(2)
   ...     client1.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
   ...                            "samp.result": {"txt": "printed sync"}})


Bind test MType for sync calls::

   >>> client1.bindReceiveCall("samp.test", test_receive_sync_call)
   >>> try:
   ...     # Sync call
   ...     print(client2.callAndWait(client1.getPublicId(),
   ...                               {"samp.mtype": "samp.test",
   ...                                "samp.params": {"txt": "Hello SYNCRO client 1!"}},
   ...                                "10"))
   ... except SAMPProxyError as e:
   ...     # If timeout expires than a SAMPProxyError is returned
   ...     print("Error ({0}): {1}".format(e.faultCode, e.faultString))
   SYNC Call: cli#2 msg#3;;cli#hub;;cli#2;;sampy::sync::call samp.test {'txt': 'Hello SYNCRO Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
   {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed sync'}}


Disconnect the clients from the hub at the end::

   >>> client1.disconnect()
   >>> client2.disconnect()
