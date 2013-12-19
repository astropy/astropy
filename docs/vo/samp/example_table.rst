.. _vo-samp-example_table:

.. doctest-skip-all

**********************************
Sending/receiving tables over SAMP
**********************************

In the following examples, we use the `TOPCAT
<http://www.star.bris.ac.uk/~mbt/topcat/>`_ tool to run a SAMP hub and to
exchange tables with. We also make use of `SAO Ds9
<http://hea-www.harvard.edu/RD/ds9>`_ since it can receive tables. Before
running the following examples, make sure that you open up TOPCAT and then Ds9.

Sending a table to TOPCAT and Ds9
=================================

The easiest way to send a VO table to TOPCAT is to make use of the
`~astropy.vo.samp.SAMPIntegratedClient` class. Once TOPCAT is open, then first
instantiate a `~astropy.vo.samp.SAMPIntegratedClient` instance and connect to
the hub::

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
``table.load.votable`` messages::

    >>> client.notifyAll(message)

The above message will actually be broadcast to all applications connected via
SAMP. For example, if we open `SAO Ds9 <http://hea-www.harvard.edu/RD/ds9>`_ in
addition to TOPCAT, and we run the above command, both applications will load
the table. We can use the ``getRegisteredClients`` method to find all the
clients connected to the hub::

    >>> client.getRegisteredClients()
    ['hub', 'c1', 'c2']
    
These IDs don't mean much, but we can find out more using::

   >>> client.getMetadata('c1')
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
this time only to TOPCAT, using::

    >>> client.notify('c1', message)

Once finished, we should make sure we disconnect from the hub:

    >>> client.disconnect()

Receiving a table from TOPCAT
=============================

To receive a table from TOPCAT, we have to set up a client that listens for
messages from the hub. As before, we instantiate a
`~astropy.vo.samp.SAMPIntegratedClient` instance and connect to the hub::

    >>> from astropy.vo.samp import SAMPIntegratedClient
    >>> client = SAMPIntegratedClient()
    >>> client.connect()

We now set up a receiver class which will handle any received message::

    >>> class Receiver(object):
    ...     def __init__(self):
    ...         self.received = False
    ...     def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
    ...         self.params = params
    ...         self.received = True

and we instantiate it:

    >>> r = Receiver()

We can now use the ``bindReceiveCall`` method to tell our receiver to listed to
all ``table.load.*`` messages::

    >>> client.bindReceiveCall("table.load.*", r.receive_call)

We can now check that the message has not been received yet::

    >>> r.received
    False
    
Let's now broadcast the table from TOPCAT. After a few seconds, we can try and check again if the message has been received::

    >>> r.received
    True
    
Success! The table URL should now be available in ``r.params['url']``, so we can do::

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
        def __init__(self):
            self.received = False
        def receive_call(self, private_key, sender_id, msg_id, mtype, params, extra):
            self.params = params
            self.received = True

    # Instantiate the receiver
    r = Receiver()

    # Listen for any instructions to load a table
    client.bindReceiveCall("table.load.*", r.receive_call)

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