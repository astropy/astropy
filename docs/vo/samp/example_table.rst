.. _vo-samp-example_table:

.. doctest-skip-all

*****************************
Exchanging tables with TOPCAT
*****************************

Sending a table to TOPCAT
=========================

The easiest way to send a VO table to TOPCAT is to make use of the
`~astropy.vo.samp.SAMPIntegratedClient` class. Once TOPCAT is open, the
following script demonstrates how to broadcast a local table over SAMP::

    from astropy.vo.samp import SAMPIntegratedClient

    # Instantiate the client and connect to the hub
    client = SAMPIntegratedClient()
    client.connect()

    # Set the parameters for the table to be sent. This should include a local
    # or remote URL, and a name for the table.
    params = {}
    params["url"] = 'file:///Users/tom/Desktop/aj285677t3_votable.xml'
    params["name"] = "Robitaille et al. (2008), Table 3"

    # Then, set the parameters for the message itself. This includes the type
    # of message (here instructing SAMP clients to load a VO table) and the
    # details of the table.
    message = {}
    message["samp.mtype"] = "table.load.votable"
    message["samp.params"] = params

    # Finally, send the message to all clients
    client.notifyAll(message)

    # Disconnect from the hub
    client.disconnect()

Note that to construct a local URL, you can also make use of ``urlparse`` as follows:

    import urlparse
    params["url"] = urlparse.urljoin('file:', os.path.abspath("aj285677t3_votable.xml"))

.. TODO: more details about sending the message only to TOPCAT

Receiving a table from TOPCAT
=============================

To receive a table from TOPCAT, we have to set up a client that listens for
messages from the hub. The following is a complete example of a simple client
that listens for requests to load a table and loads it when received::

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

Sending a table from TOPCAT then results in the client printing out e.g.

    Downloading http://127.0.0.1:2525/dynamic/4/t12.vot [Done]
               col1             col2     col3    col4     col5    col6 col7  col8 col9 col10
    ------------------------- -------- ------- -------- -------- ----- ---- ----- ---- -----
    SSTGLMC G000.0046+01.1431   0.0046  1.1432 265.2992 -28.3321  6.67 5.04  6.89 5.22     N
    SSTGLMC G000.0106-00.7315   0.0106 -0.7314 267.1274 -29.3063  7.18 6.07   nan 5.17     Y
    SSTGLMC G000.0110-01.0237   0.0110 -1.0236 267.4151 -29.4564  8.32 6.30  8.34 6.32     N
    ...
