Add ``make_rgb()``, ``make_log_rgb()``, and ``make_linear_rgb()``, convenience 
functions for creating RGB images with independent scaling on each filter. 
Refactors ``make_lupton_rgb()`` to work with instances of subclasses of 
``BaseStretch``, including the new Lupton-specific classes 
``AsinhLuptonStretch`` and ``AsinhZscaleLuptonStretch``.