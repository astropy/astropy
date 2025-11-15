Fixed misleading ``AttributeError`` messages when accessing properties in
subclassed ``SkyCoord`` objects. Previously, if a property raised an
``AttributeError`` for a missing attribute, the error message would
incorrectly state that the property itself didn't exist, rather than the
attribute that was missing within the property.
