# Licensed under a 3-clause BSD style license - see LICENSE.rst
#This module implements the base NDData class.

__all__ = ['NDData']

import numpy as np

from .masks import NDMask, BoolMask
from .errors import NDError, SDError

from astropy.config import logger

def merge_meta(meta1, meta2):
    #Merging meta information and removing all duplicate keys -- warning what keys were removed
    #should be in NDData somewhere
    meta1_keys = meta1.viewkeys()
    meta2_keys = meta2.viewkeys()
    
    
    duplicates = meta1_keys & meta2_keys
    
    if len(duplicates) > 0:
        logger.warn('Removing duplicate keys found in meta data: ' + ','.join(duplicates))
                    
    new_meta = copy.deepcopy(meta1)
    new_meta.update(copy.deepcopy(meta2))
    for key in duplicates:
        del new_meta[key]
    
    return new_meta



def nddata_operation(func):
    #used as a decorator for the arithmetic of spectra
    def convert_operands(self, operand):
        
        #checking if they have the same wcs and units
        if isinstance(operand, self.__class__):
            #check if wcs is the same ! Can't be implemented yet due to missing wcs (in spectrum1d it checks the dispersion)
            
            data = operand.data
            mask = operand.mask
            error = operand.error
            meta = operand.meta
        #for a scalar the data is the scalar and the error and mask and meta are None/ {}
        elif np.isscalar(operand):
            data = operand
            mask = None
            error = None
            meta = {}
            
        else:
            raise ValueError("unsupported operand type(s) for operation: %s and %s" %
                             (type(self), type(operand)))
        
        return func(self, data, error, mask, meta)
        
    return convert_operands


class NDData(object):
    """Superclass for Astropy data.

    `NDData` provides a superclass for all array-based data. The key
    distinction from raw numpy arrays is the presence of additional metadata
    such as error arrays, bad pixel masks, units and/or a coordinate system.

    Parameters
    -----------
    data : `~numpy.ndarray`
        The actual data contained in this `NDData` object.
        
    error : `~numpy.ndarray`, optional
        Error of the data. This should be interpreted as a 1-sigma error (e.g,
        square root of the variance), under the assumption of Gaussian errors.
        Must be a shape that can be broadcast onto `data`.

         .. warning::
             The physical interpretation of the `error` array may change in the
             future, as it has not been intensively discussed. For now assume
             the above description holds, using an `error` property if
             necessary, but feel free to use the most convinient internal
             representation in subclasses

    mask : `~numpy.ndarray`, optional
        Masking of the data; Should be False/0 (or the empty string) where the
        data is *valid*.  All other values indicate that the value should be
        masked. Must be a shape that can be broadcast onto `data`.
        
    wcs : undefined, optional
        WCS-object containing the world coordinate system for the data.

        .. warning::
            This is not yet defind because the discussion of how best to
            represent this class's WCS system generically is still under
            consideration. For now just leave it as None

    meta : `dict`-like object, optional
        Metadata for this object.  "Metadata" here means all information that
        is included with this object but not part of any other attribute 
        of this particular object.  e.g., creation date, unique identifier, 
        simulation parameters, exposure time, telescope name, etc.
        
    units : undefined, optional
        The units of the data.

        .. warning::
            The units scheme is under development. For now, just supply a
            string when relevant - the units system will likely be compatible
            with providing strings to initialize itself.

    copy : bool, optional
        If True, the array will be *copied* from the provided `data`, otherwise
        it will be referenced if possible (see `numpy.array` :attr:`copy`
        argument for details).
        
    validate : bool, optional
        If False, no type or shape-checking or array conversion will occur.
        Note that if `validate` is False, :attr:`copy` will be ignored.

    Raises
    ------
    ValueError
        If the `error` or `mask` inputs cannot be broadcast (e.g., match
        shape) onto `data`.

    """
    def __init__(self, data, error=None, mask=None, wcs=None, meta=None,
                 units=None, copy=True, validate=True):
        if validate:
            self.data = np.array(data, subok=True, copy=copy)
            self.error = error
            self.mask = mask

            self._validate_mask_and_error()

            self.wcs = wcs
            self.units = units
            if meta is None:
                self.meta = {}
            else:
                self.meta = dict(meta)  # makes a *copy* of the passed-in meta
        else:
            self.data = np.array(data, copy=copy)
            self.error = np.array(error, copy=copy)
            self.mask = np.array(mask, copy=copy)
            self.meta = copy.deepcopy(meta) if copy else meta
            
            #what about copy for units and wcs?
            self.units = units
            self.wcs = wcs
    def _validate_mask_and_error(self):
        """
        Raises ValueError if they don't match (using ~numpy.broadcast)
        """
        if self.mask is not None:
            if isinstance(self.mask, np.ndarray):
                self.mask = nddata.BoolMask(mask)
            elif isinstance(self.mask, nddata.NDMask):
                pass
            else:
                raise ValueError('Mask not of appropriate type')
            
        try:
            if self.mask is not None:
                np.broadcast(self.data, self.mask)
            maskmatch = True
        except ValueError:
            maskmatch = False

        if self.error is not None:
            if isinstance(self.error, np.ndarray):
                    self.error = nddata.SDError(self.error)
            elif isinstance(self.error, nddata.NDError):
                pass
            else:
                raise ValueError('Error not of appropriate type')

        try:
            if self.error is not None:
                np.broadcast(self.data, self.error)
            errmatch = True
        except ValueError:
            errmatch = False

        if not errmatch and not maskmatch:
            raise ValueError('NDData error and mask do not match data')
        elif not errmatch:
            raise ValueError('NDData error does not match data')
        elif not maskmatch:
            raise ValueError('NDData mask does not match data')

    @property
    def boolmask(self):
        """
        The mask as a boolean array (or None if the mask is None).

        This mask is True where the data is *valid*, and False where the data
        should be *masked*.  This is the opposite of the convention used for
        `mask`, but allows simple retrieval of the unmasked data points as
        ``ndd.data[ndd.boolmask]``.
        """
        if self.mask is None:
            return None
        else:
            dtchar = self.mask.dtype.char
            if dtchar is 'U' or dtchar is 'S':
                return self.mask == ''
            else:
                return ~self.mask.astype(bool)

    @property
    def shape(self):
        """
        shape tuple of this object's data.
        """
        return self.data.shape

    @property
    def size(self):
        """
        integer size of this object's data.
        """
        return self.data.size

    @property
    def dtype(self):
        """
        `numpy.dtype` of this object's data.
        """
        return self.data.dtype

    @property
    def ndim(self):
        """
        integer dimensions of this object's data
        """
        return self.data.ndim
        
    @nddata_operation 
    def __add__(self, operand_data, operand_error, operand_mask, operand_meta):
        new_data = self.data + operand_data
        
        if self.error is None and operand_error is None:
            new_error = None
        elif self.error is None:
            new_error = operand_error.copy()
        else:
            new_error = self.error.error_add(self.data, operand_error, operand_data, new_data)
        
        if self.mask is None and operand_mask is None:
            new_mask = None
        elif self.mask is None:
            new_mask = operand_mask.copy()
        else:
            new_mask = self.mask.mask_add(operand_mask)
        
        new_meta = merge_meta(self.meta, operand_meta)
        
        return self.__class__(new_data,
                              #!!! What if it's a WCS
                              self.dispersion.copy(),
                              error=new_error,
                              mask=new_mask,
                              meta=new_meta)
