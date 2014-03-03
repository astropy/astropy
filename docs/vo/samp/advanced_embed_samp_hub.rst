.. include:: references.txt

.. doctest-skip-all

Embedding a SAMP hub in a GUI
-----------------------------

Overview
^^^^^^^^

If you wish to embed a SAMP hub in your Python GUI tool, you will need to start
the hub programmatically using::

    from astropy.vo.samp import SAMPHubServer
    hub = SAMPHubServer()
    hub.start()

This launches the hub in a thread and is non-blocking. If you are not
interested in connections from web SAMP clients, then you can simply use::

    from astropy.vo.samp import SAMPHubServer
    hub = SAMPHubServer(web_profile=False)
    hub.start()

and this should be all you need to do. However, if you want to keep the web
profile active, there is an additional consideration, which is that when a web
SAMP client connects, you will need to ask the user whether they accept the
connection (for security reasons). By default, the confirmation message is a
text-based message in the terminal, but if you have a GUI tool, you will
instead likely want to open a GUI dialog.

To do this, you will need to define a class that handles the dialog,
and you should then pass an **instance** of the class to
|SAMPHubServer| (not the class itself).  This class should inherit
from `astropy.vo.samp.WebProfileDialog` and add the following:

    1) It should have a GUI timer callback that periodically calls
       ``WebProfileDialog.handle_queue`` (available as
       ``self.handle_queue``).

    2) Implement a ``show_dialog`` method to display a consent dialog.
       It should take the following arguments:

           - ``samp_name``: The name of the application making the request.

           - ``details``: A dictionary of details about the client
             making the request. The only key in this dictionary required by
             the SAMP standard is ``samp.name`` which gives the name of the
             client making the request.

           - ``client``: A hostname, port pair containing the client
             address.

           - ``origin``: A string containing the origin of the
             request.

    3) Based on the user response, the ``show_dialog`` should call
       ``WebProfileDialog.consent`` or ``WebProfileDialog.reject``.
       This may, in some cases, be the result of another GUI callback.

Example of embedding a SAMP hub in a Tk application
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following code is a full example of a simple Tk application that watches
for web SAMP connections and opens the appropriate dialog::


    import Tkinter as tk
    import tkMessageBox

    from astropy.vo.samp import SAMPHubServer
    from astropy.vo.samp.hub import WebProfileDialog

    MESSAGE = """
    A Web application which declares to be

    Name: {name}
    Origin: {origin}

    is requesting to be registered with the SAMP Hub.  Pay attention
    that if you permit its registration, such application will acquire
    all current user privileges, like file read/write.

    Do you give your consent?
    """

    class TkWebProfileDialog(WebProfileDialog):
        def __init__(self, root):
            self.root = root
            self.wait_for_dialog()

        def wait_for_dialog(self):
            self.handle_queue()
            self.root.after(100, self.wait_for_dialog)

        def show_dialog(self, samp_name, details, client, origin):
            text = MESSAGE.format(name=samp_name, origin=origin)

            response = tkMessageBox.askyesno(
                'SAMP Hub', text,
                default=tkMessageBox.NO)

            if response:
                self.consent()
            else:
                self.reject()

    # Start up Tk application
    root = tk.Tk()
    tk.Label(root, text="Example SAMP Tk application",
             font=("Helvetica", 36), justify=tk.CENTER).pack(pady=200)
    root.geometry("500x500")
    root.update()

    # Start up SAMP hub
    h = SAMPHubServer(web_profile_dialog=TkWebProfileDialog(root))
    h.start()

    try:
        # Main GUI loop
        root.mainloop()
    except KeyboardInterrupt:
        pass

    h.stop()

If you run the above script, a window will open saying "Example SAMP Tk
application". If you then go to the following page for example:

http://astrojs.github.io/sampjs/examples/pinger.html

and click on the Ping button, you will see the dialog open in the Tk
application. Once you click on 'CONFIRM', future 'Ping' calls will no longer
bring up the dialog.
