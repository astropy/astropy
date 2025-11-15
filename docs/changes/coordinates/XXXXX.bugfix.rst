Fixed misleading ``AttributeError`` messages when accessing properties in subclassed
``SkyCoord`` objects. Previously, if a custom property tried to access a non-existent
attribute, the error message would incorrectly claim that the property itself didn't exist.
Now the error message clearly indicates that the property exists but raised an ``AttributeError``,
helping developers identify the actual missing attribute.
