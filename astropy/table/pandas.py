ascii_coded = 'Ò♙♙♙♙♙♙♙♙♌♐♐♌♙♙♙♙♙♙♌♌♙♙Ò♙♙♙♙♙♙♙♘♐♐♐♈♙♙♙♙♙♌♐♐♐♔Ò♙♙♌♈♙♙♌♐♈♈♙♙♙♙♙♙♙♙♈♐♐♙Ò♙♐♙♙♙♐♐♙♙♙♙♙♙♙♙♙♙♙♙♙♙♙Ò♐♔♙♙♘♐♐♙♙♌♐♐♔♙♙♌♌♌♙♙♙♌Ò♐♐♙♙♘♐♐♌♙♈♐♈♙♙♙♈♐♐♙♙♘♔Ò♐♐♌♙♘♐♐♐♌♌♙♙♌♌♌♙♈♈♙♌♐♐Ò♘♐♐♐♌♐♐♐♐♐♐♌♙♈♙♌♐♐♐♐♐♔Ò♘♐♐♐♐♐♐♐♐♐♐♐♐♈♈♐♐♐♐♐♐♙Ò♙♘♐♐♐♐♈♐♐♐♐♐♐♙♙♐♐♐♐♐♙♙Ò♙♙♙♈♈♈♙♙♐♐♐♐♐♔♙♐♐♐♐♈♙♙Ò♙♙♙♙♙♙♙♙♙♈♈♐♐♐♙♈♈♈♙♙♙♙Ò'
ascii_uncoded = ''.join([chr(ord(c)-200) for c in ascii_coded])

try:
    from IPython import display
    im = display.Image(url='https://media.giphy.com/media/e24Q8FKE2mxRS/giphy.gif')
    display.display(im)
except ImportError:
    pass

print(ascii_uncoded)
