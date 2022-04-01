# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import tkinter
from tkinter import ttk

def make_modal(text='CALCULATIONS\nCORRECT',
               title='Are you plugged into Pleiades?',
               blink_ms=750,
               geometry='550x155',
               font=('Helvetica', 48)):
    tkroot = tkinter.Tk()
    tkroot.title(title)
    tkroot.geometry(geometry)
    tkroot.configure(bg='black')

    s = ttk.Style()
    s.configure('my.TButton', font=font,
                              background='black',
                              foreground='white',
                              justify='center',
                              anchor='center')

    frm = ttk.Frame(tkroot, padding=5)
    frm.pack(side=tkinter.LEFT, expand=True, fill='both')

    last_afterjob = [-1]
    def destroy():
        if last_afterjob[0] != -1:
            tkroot.after_cancel(last_afterjob[0])
        tkroot.destroy()
    btn = ttk.Button(frm, text='', command=destroy, style='my.TButton')

    btn.pack(expand=True, fill='both')

    def blink_text():
        if btn['text'] == text:
            btn['text'] = ''
        else:
            btn['text'] = text
        last_afterjob[0] = tkroot.after(blink_ms, blink_text)
    blink_text()

    tkroot.wait_visibility()
    tkroot.grab_set_global()

    tkroot.mainloop()