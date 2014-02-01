.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example_clients:


Communication between integrated clients objects
------------------------------------------------

As shown in :doc:`example_table_image`, the |SAMPIntegratedClient| class can be
used to communicate with other SAMP-enabled tools such as `TOPCAT
<http://www.star.bris.ac.uk/~mbt/topcat/>`_, `SAO Ds9
<http://hea-www.harvard.edu/RD/ds9>`_, or `Aladin Desktop
<http://aladin.u-strasbg.fr>`_.

In this section, we look at how we can set up two |SAMPIntegratedClient|
instances and communicate between them.

First, start up a SAMP hub as described in :doc:`example_hub`.

Next, we create two clients and connect them to the hub::

   >>> client1 = samp.SAMPIntegratedClient(name="Client 1", description="Test Client 1",
   ...                                     metadata = {"client1.version":"0.01"})
   >>> client2 = samp.SAMPIntegratedClient(name="Client 2", description="Test Client 2",
   ...                                     metadata = {"client2.version":"0.25"})
   >>> client1.connect()
   >>> client2.connect()

We now define functions to call when receiving a notification, call or response::

   >>> def test_receive_notification(private_key, sender_id, mtype, params, extra):
   ...     print("Notification:", private_key, sender_id, mtype, params, extra)

   >>> def test_receive_call(private_key, sender_id, msg_id, mtype, params, extra):
   ...     print("Call:", private_key, sender_id, msg_id, mtype, params, extra)
   ...     client1.ereply(msg_id, SAMP_STATUS_OK, result = {"txt": "printed"})

   >>> def test_receive_response(private_key, sender_id, msg_id, response):
   ...     print("Response:", private_key, sender_id, msg_id, response)

We subscribe client 1 to ``"samp.*"`` and ``"samp.app.*"`` and bind them to the
related functions::

   >>> client1.bind_receive_notification("samp.app.*", test_receive_notification)
   >>> client1.bind_receive_call("samp.app.*", test_receive_call)

We now bind message tags received by client 2 to suitable functions::

   >>> client2.bind_receive_response("my-dummy-print", test_receive_response)
   >>> client2.bind_receive_response("my-dummy-print-specific", test_receive_response)

We are now ready to test out the clients and callback functions. Client 2
notifies all clients using the "samp.app.echo" message type via the hub::

   >>> client2.enotify_all("samp.app.echo", txt="Hello world!")
   ['cli#2']
   Notification: 0d7f4500225981c104a197c7666a8e4e cli#2 samp.app.echo {'txt':
   'Hello world!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}

We can also find a dictionary giving the clients that would currently receive
``samp.app.echo`` messages::

   >>> print(client2.getSubscribedClients("samp.app.echo"))
   {'cli#2': {}}

Client 2 calls all clients with the ``"samp.app.echo"`` message type using
``"my-dummy-print"`` as a message-tag::

   >>> print(client2.call_all("my-dummy-print",
   ...                        {"samp.mtype": "samp.app.echo",
   ...                         "samp.params": {"txt": "Hello world!"}}))
   {'cli#1': 'msg#1;;cli#hub;;cli#2;;my-dummy-print'}
   Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2
   msg#1;;cli#hub;;cli#2;;my-dummy-print samp.app.echo {'txt': 'Hello world!'}
   {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
   Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print
   {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}

Client 2 then calls client 1 using the ``"samp.app.echo"`` message type,
tagging the message as ``"my-dummy-print-specific"``::

   >>> try:
   ...     print(client2.call(client1.getPublicId(),
   ...                        "my-dummy-print-specific",
   ...                        {"samp.mtype": "samp.app.echo",
   ...                         "samp.params": {"txt": "Hello client 1!"}}))
   ... except SAMPProxyError as e:
   ...     print("Error ({0}): {1}".format(e.faultCode, e.faultString))
   msg#2;;cli#hub;;cli#2;;my-dummy-print-specific
   Call: 8c8eb53178cb95e168ab17ec4eac2353 cli#2
   msg#2;;cli#hub;;cli#2;;my-dummy-print-specific samp.app.echo {'txt': 'Hello
   Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user': 'unknown'}
   Response: d0a28636321948ccff45edaf40888c54 cli#1 my-dummy-print-specific
   {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed'}}

We can now define a function called to test synchronous calls::

   >>> def test_receive_sync_call(private_key, sender_id, msg_id, mtype, params, extra):
   ...     import time
   ...     print("SYNC Call:", sender_id, msg_id, mtype, params, extra)
   ...     time.sleep(2)
   ...     client1.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
   ...                            "samp.result": {"txt": "printed sync"}})

We now bind the ``samp.test`` message type to ``test_receive_sync_call``::

   >>> client1.bind_receive_call("samp.test", test_receive_sync_call)
   >>> try:
   ...     # Sync call
   ...     print(client2.call_and_wait(client1.getPublicId(),
   ...                                 {"samp.mtype": "samp.test",
   ...                                  "samp.params": {"txt": "Hello SYNCRO client 1!"}},
   ...                                  "10"))
   ... except SAMPProxyError as e:
   ...     # If timeout expires than a SAMPProxyError is returned
   ...     print("Error ({0}): {1}".format(e.faultCode, e.faultString))
   SYNC Call: cli#2 msg#3;;cli#hub;;cli#2;;sampy::sync::call samp.test {'txt':
   'Hello SYNCRO Cli 1!'} {'host': 'antigone.lambrate.inaf.it', 'user':
   'unknown'}
   {'samp.status': 'samp.ok', 'samp.result': {'txt': 'printed sync'}}

Finally, we disconnect the clients from the hub at the end::

   >>> client1.disconnect()
   >>> client2.disconnect()
