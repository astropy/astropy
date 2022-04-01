# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['make_modal']

try:
    import tkinter
    from tkinter import ttk
except ImportError:
    tkinter = None

if tkinter is None:
    def make_modal(text, title, blink_ms, geometry, font, center_window=True):
        print(text)
else:
    def make_modal(text, title, blink_ms, geometry, font, center_window=True):
        tkroot = tkinter.Tk()

        if center_window:
            screen_width = tkroot.winfo_screenwidth()
            screen_height = tkroot.winfo_screenheight()
            width, height = geometry.split('x')
            x = int(screen_width/2 - int(width)/2)
            y = int(screen_height/2 - int(height)/2)
            tkroot.geometry(geometry + f'+{x}+{y}')
        else:
            tkroot.geometry(geometry)

        tkroot.title(title)
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