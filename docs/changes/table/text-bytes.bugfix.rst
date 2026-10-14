Avoid implicit NumPy bytes/text promotion when importing masked string columns
from DataFrames, assigning mixed bytes and text to byte columns, and joining or
stacking byte and Unicode columns. Mixed table columns now use UTF-8 decoding.
