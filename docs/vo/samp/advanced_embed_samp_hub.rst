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
SAMP client connects, you will need to ask the user whehter they accept the
connection (for security reasons). By default, the confirmation messge is a
text-based message in the terminal, but if you have a GUI tool, you will
instead likely want to open a GUI dialog.

To do this, there are two things you need to do. First, disable the default text-based dialog::

   hub.set_web_profile_text_dialog(False)

and secondly, run the following command to get the queues which will allow your
GUI to communicate with the SAMP hub::

    queue_request, queue_result = hub.get_web_profile_dialog_queues()

These are both standard instances of `Queue.Queue`. If a web client asks to
connect to the SAMP hub, the ``queue_request`` queue is passed a tuple
containing: a dictionary with details of the client that is connecting, a
string containing the client address, and a string containing the origin of the
request. The following is an example of such a tuple::

    ({'samp.name': 'Monitor'}, ['127.0.0.1', 58105], 'http://astrojs.github.io')

Once your GUI sees such an object in the queue, it should open a dialog asking
the user whether to approve the connection. If the user accepts, you should
pass `True` to ``queue_result``, otherwise you should pass `False`.

Example of embedding a SAMP hub in a Tk application
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following code is a full example of a simple Tk application that watches
for web SAMP connections and opens the appropriate dialog::

    import Queue
    import Tkinter as tk

    from astropy.vo.samp import SAMPHubServer

    MESSAGE = """
    A Web application which declares to be

    Name: {name}
    Origin: {origin}

    is requesting to be registered with the SAMP Hub.
    Pay attention that if you permit its registration, such
    application will acquire all current user privileges, like
    file read/write.

    Do you give your consent?
    """


    class WebProfileDialog(object):

        def __init__(self, root, queue_in, queue_out):

            self.root = root

            self._queue_in = queue_in
            self._queue_out = queue_out

            self.wait_for_dialog()

        def wait_for_dialog(self):
            try:
                request = self._queue_in.get_nowait()
            except Queue.Empty:
                pass
            else:
                self.show_dialog(request)
            self.root.after(100, self.wait_for_dialog)

        def show_dialog(self, request):

            self.window = tk.Toplevel(self.root)

            self.window.title("SAMP Hub")

            self.label = tk.Label(self.window, font=("Helvetica", 14),
                                  fg="red", justify=tk.CENTER)
            self.label.pack(padx=5, pady=5, expand=1, fill=tk.X,)

            a = tk.Button(self.window, text="CONSENT", command=self.consent)
            a.pack(padx=5, pady=5, expand=1, fill=tk.X, side=tk.LEFT)

            r = tk.Button(self.window, text="REJECT", command=self.reject)
            r.pack(padx=5, pady=5, expand=1, fill=tk.X, side=tk.RIGHT)

            self.window.protocol("WM_DELETE_WINDOW", self.reject)

            if isinstance(request[0], str):  # To support the old protocol version
                samp_name = request[0]
            else:
                samp_name = request[0]["samp.name"]

            text = MESSAGE.format(name=samp_name, origin=request[2])

            self.label.configure(text=text)
            self.window.update()

        def consent(self):
            self._queue_out.put(True)
            self.window.destroy()

        def reject(self):
            self._queue_out.put(False)
            self.window.destroy()

    # Start up Tk application
    root = tk.Tk()
    tk.Label(root, text="Example SAMP Tk application",
             font=("Helvetica", 36), justify=tk.CENTER).pack(pady=200)
    root.geometry("500x500")
    root.update()

    # Start up SAMP hub
    h = SAMPHubServer()
    h.set_web_profile_text_dialog(False)
    queue_in, queue_out = h.get_web_profile_dialog_queues()
    h.start()

    # Prepare dialog in case it is needed
    d = WebProfileDialog(root, queue_in, queue_out)

    # Main GUI loop
    root.mainloop()

If you run the above script, a window will open saying "Example SAMP Tk
application". If you then go to the following page for example:

http://astrojs.github.io/sampjs/examples/pinger.html

and click on the Ping button, you will see the dialog open in the Tk
application. Once you click on 'CONFIRM', future 'Ping' calls will no longer
bring up the dialog.
