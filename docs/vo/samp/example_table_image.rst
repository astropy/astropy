.. include:: references.txt

.. doctest-skip-all

.. _vo-samp-example-table-image:

Sending/receiving tables and images over SAMP
---------------------------------------------

In the following examples, we make use of:

* `TOPCAT <http://www.star.bris.ac.uk/~mbt/topcat/>`_, which is a tool to
  explore tabular data.
* `SAO Ds9 <http://hea-www.harvard.edu/RD/ds9>`_, which is an image
  visualization tool, which can also overplot catalogs.
* `Aladin Desktop <http://aladin.u-strasbg.fr>`_, which is another tool that
  can visualize images and catalogs.

TOPCAT and Aladin will run a SAMP Hub is none is found, so for the following
examples you can either start up one of these applications first, or you can
start up the `astropy.vo.samp` hub. You can start this using the following
command::

    $ samp_hub

Sending a table to TOPCAT and Ds9
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to send a VO table to TOPCAT is to make use of the
|SAMPIntegratedClient| class. Once TOPCAT is open, then first instantiate a
|SAMPIntegratedClient| instance and connect to the hub::

    >>> from astropy.vo.samp import SAMPIntegratedClient
    >>> client = SAMPIntegratedClient()
    >>> client.connect()

Next, we have to set up a dictionary that contains details about the table to
send. This should include ``url``, which is the URL to the file, and ``name``,
which is a human-readable name for the table. The URL can be a local URL
(starting with ``file:///``)::

    >>> params = {}
    >>> params["url"] = 'file:///Users/tom/Desktop/aj285677t3_votable.xml'
    >>> params["name"] = "Robitaille et al. (2008), Table 3"

.. note:: To construct a local URL, you can also make use of ``urlparse`` as
          follows::

                >>> import urlparse
                >>> params["url"] = urlparse.urljoin('file:', os.path.abspath("aj285677t3_votable.xml"))

Now we can set up the message itself. This includes the type of message (here
we use ``table.load.votable`` which indicates that a VO table should be loaded,
and the details of the table that we set above::

    >>> message = {}
    >>> message["samp.mtype"] = "table.load.votable"
    >>> message["samp.params"] = params

Finally, we can broadcast this to all clients that are listening for
``table.load.votable`` messages using
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.notify_all`::

    >>> client.notify_all(message)

The above message will actually be broadcast to all applications connected via
SAMP. For example, if we open `SAO Ds9 <http://hea-www.harvard.edu/RD/ds9>`_ in
addition to TOPCAT, and we run the above command, both applications will load
the table. We can use the
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.get_registered_clients` method to
find all the clients connected to the hub::

    >>> client.get_registered_clients()
    ['hub', 'c1', 'c2']

These IDs don't mean much, but we can find out more using::

   >>> client.get_metadata('c1')
   {'author.affiliation': 'Astrophysics Group, Bristol University',
    'author.email': 'm.b.taylor@bristol.ac.uk',
    'author.name': 'Mark Taylor',
    'home.page': 'http://www.starlink.ac.uk/topcat/',
    'samp.description.text': 'Tool for OPerations on Catalogues And Tables',
    'samp.documentation.url': 'http://127.0.0.1:2525/doc/sun253/index.html',
    'samp.icon.url': 'http://127.0.0.1:2525/doc/images/tc_sok.gif',
    'samp.name': 'topcat',
    'topcat.version': '4.0-1'}

We can see that ``c1`` is the TOPCAT client. We can now re-send the data, but
this time only to TOPCAT, using the
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.notify` method::

    >>> client.notify('c1', message)

Once finished, we should make sure we disconnect from the hub::

    >>> client.disconnect()

Receiving a table from TOPCAT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To receive a table from TOPCAT, we have to set up a client that listens for
messages from the hub. As before, we instantiate a |SAMPIntegratedClient|
instance and connect to the hub::

    >>> from astropy.vo.samp import SAMPIntegratedClient
    >>> client = SAMPIntegratedClient()
    >>> client.connect()

We now set up a receiver class which will handle any received message. We need
to take care to write handlers for both notifications and calls (the difference
between the two being that calls expect a reply)::

    >>> class Receiver(object):
    ...     def __init__(self, client):
    ...         self.client = client
    ...         self.received = False
    ...     def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
    ...         self.params = params
    ...         self.received = True
    ...         self.client.reply(msg_id, {"samp.status": "samp.ok", "samp.result": {}})
    ...     def receive_notification(self, private_key, sender_id, mtype, params, extra):
    ...         self.params = params
    ...         self.received = True

and we instantiate it:

    >>> r = Receiver(client)

We can now use the
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.bind_receive_call` and
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.bind_receive_notification` methods
to tell our receiver to listen to all ``table.load.votable`` messages::

    >>> client.bind_receive_call("table.load.votable", r.receive_call)
    >>> client.bind_receive_notification("table.load.votable", r.receive_notification)

We can now check that the message has not been received yet::

    >>> r.received
    False

Let's now broadcast the table from TOPCAT. After a few seconds, we can try and
check again if the message has been received::

    >>> r.received
    True

Success! The table URL should now be available in ``r.params['url']``, so we can do::

    >>> from astropy.table import Table
    >>> t = Table.read(r.params['url'])
    Downloading http://127.0.0.1:2525/dynamic/4/t12.vot [Done]
    >>> t
               col1             col2     col3    col4     col5    col6 col7  col8 col9 col10
    ------------------------- -------- ------- -------- -------- ----- ---- ----- ---- -----
    SSTGLMC G000.0046+01.1431   0.0046  1.1432 265.2992 -28.3321  6.67 5.04  6.89 5.22     N
    SSTGLMC G000.0106-00.7315   0.0106 -0.7314 267.1274 -29.3063  7.18 6.07   nan 5.17     Y
    SSTGLMC G000.0110-01.0237   0.0110 -1.0236 267.4151 -29.4564  8.32 6.30  8.34 6.32     N
    ...

As before, we should remember to disconnect from the hub once we are done::

    >>> client.disconnect()

The following is a full example of a script that can be used to receive and
read a table. It includes a loop that waits until the message is received, and
reads the table once it has::

    import time

    from astropy.vo.samp import SAMPIntegratedClient
    from astropy.table import Table

     # Instantiate the client and connect to the hub
    client=SAMPIntegratedClient()
    client.connect()

    # Set up a receiver class
    class Receiver(object):
        def __init__(self, client):
            self.client = client
            self.received = False
        def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
            self.params = params
            self.received = True
            self.client.reply(msg_id, {"samp.status": "samp.ok", "samp.result": {}})
        def receive_notification(self, private_key, sender_id, mtype, params, extra):
            self.params = params
            self.received = True

    # Instantiate the receiver
    r = Receiver(client)

    # Listen for any instructions to load a table
    client.bind_receive_call("table.load.votable", r.receive_call)
    client.bind_receive_notification("table.load.votable", r.receive_notification)

    # We now run the loop to wait for the message in a try/finally block so that if
    # the program is interrupted e.g. by control-C, the client terminates
    # gracefully.

    try:

        # We test every 0.1s to see if the hub has sent a message
        while True:
            time.sleep(0.1)
            if r.received:
                t = Table.read(r.params['url'])
                break

    finally:

        client.disconnect()

    # Print out table
    print t

Sending an image to Ds9 and Aladin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As for tables, the easiest way to send a FITS image over SAMP is to make use of
the |SAMPIntegratedClient| class. Once Aladin or Ds9 are open, then first
instantiate a |SAMPIntegratedClient| instance and connect to the hub as before::

    >>> from astropy.vo.samp import SAMPIntegratedClient
    >>> client = SAMPIntegratedClient()
    >>> client.connect()

Next, we have to set up a dictionary that contains details about the image to
send. This should include ``url``, which is the URL to the file, and ``name``,
which is a human-readable name for the table. The URL can be a local URL
(starting with ``file:///``)::

    >>> params = {}
    >>> params["url"] = 'file:///Users/tom/Desktop/MSX_E.fits'
    >>> params["name"] = "MSX Band E Image of the Galactic Center"

See `Sending a table to TOPCAT and Ds9`_ for an example of how to construct local URLs
more easily. Now we can set up the message itself. This includes the type of
message (here we use ``image.load.fits`` which indicates that a FITS image
should be loaded, and the details of the table that we set above::

    >>> message = {}
    >>> message["samp.mtype"] = "image.load.fits"
    >>> message["samp.params"] = params

Finally, we can broadcast this to all clients that are listening for
``table.load.votable`` messages::

    >>> client.notify_all(message)

As for `Sending a table to TOPCAT and Ds9`_, the
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.notify_all`
method will broadcast the image to all listening clients, and as for tables it
is possible to instead use the
:meth:`~astropy.vo.samp.integrated_client.SAMPIntegratedClient.notify` method
to send it to a specific client.

Once finished, we should make sure we disconnect from the hub::

    >>> client.disconnect()

Receiving a table from Ds9 or Aladin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Receiving images over SAMP is identical to `Receiving a table from TOPCAT`_,
with the execption that the message type should be ``image.load.fits`` instead
of ``table.load.votable``. Once the URL has been received, the FITS image can
be opened with::

    >>> from astropy.io import fits
    >>> fits.open(r.params['url'])

