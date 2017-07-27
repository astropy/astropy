ascii_coded = 'Ò♙♙♙♙♙♙♙♙♌♐♐♌♙♙♙♙♙♙♌♌♙♙Ò♙♙♙♙♙♙♙♘♐♐♐♈♙♙♙♙♙♌♐♐♐♔Ò♙♙♌♈♙♙♌♐♈♈♙♙♙♙♙♙♙♙♈♐♐♙Ò♙♐♙♙♙♐♐♙♙♙♙♙♙♙♙♙♙♙♙♙♙♙Ò♐♔♙♙♘♐♐♙♙♌♐♐♔♙♙♌♌♌♙♙♙♌Ò♐♐♙♙♘♐♐♌♙♈♐♈♙♙♙♈♐♐♙♙♘♔Ò♐♐♌♙♘♐♐♐♌♌♙♙♌♌♌♙♈♈♙♌♐♐Ò♘♐♐♐♌♐♐♐♐♐♐♌♙♈♙♌♐♐♐♐♐♔Ò♘♐♐♐♐♐♐♐♐♐♐♐♐♈♈♐♐♐♐♐♐♙Ò♙♘♐♐♐♐♈♐♐♐♐♐♐♙♙♐♐♐♐♐♙♙Ò♙♙♙♈♈♈♙♙♐♐♐♐♐♔♙♐♐♐♐♈♙♙Ò♙♙♙♙♙♙♙♙♙♈♈♐♐♐♙♈♈♈♙♙♙♙Ò'
ascii_uncoded = ''.join([chr(ord(c)-200) for c in ascii_coded])
url = 'https://media.giphy.com/media/e24Q8FKE2mxRS/giphy.gif'

try:
    from IPython import display

    class ImageWithBackup(display.Image):
        def __init__(self, url, backup_text):
            super().__init__(url=url)
            self.backup_text = backup_text
        def __repr__(self):
            if self.backup_text is None:
                return super().__repr__()
            else:
                return self.backup_text


    im = ImageWithBackup(url, ascii_uncoded)
    display.display(im)
except ImportError:
    print(ascii_uncoded)
