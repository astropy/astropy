import abc

__all__ = ['BaseLowLevelWCS']


class BaseLowLevelWCS(metaclass=abc.ABCMeta):
    """
    Abstract base class for the low-level WCS interface described in APE 14
    (https://doi.org/10.5281/zenodo.1188875)
    """

    @abc.abstractproperty
    def pixel_n_dim(self):
        """
        The number of axes in the pixel coordinate system
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def world_n_dim(self):
        """
        The number of axes in the world coordinate system
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def pixel_shape(self):
        """
        The shape of the data that the WCS applies to as a tuple of
        length ``pixel_n_dim`` (optional).

        If the WCS is valid in the context of a dataset with a particular
        shape, then this property can be used to store the shape of the
        data. This can be used for example if implementing slicing of WCS
        objects. This is an optional property, and it should return `None`
        if a shape is not known or relevant.
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def pixel_bounds(self):
        """
        The bounds (in pixel coordinates) inside which the WCS is defined,
        as a list with ``pixel_n_dim`` ``(min, max)`` tuples (optional).

        WCS solutions are sometimes only guaranteed to be accurate within a
        certain range of pixel values, for example when definining a WCS
        that includes fitted distortions. This is an optional property, and
        it should return `None` if a shape is not known or relevant.
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def world_axis_physical_types(self):
        """
        Returns an iterable of strings describing the physical type for each
        world axis. They should be names from the VO UCD1+ controlled
        Vocabulary (http://www.ivoa.net/documents/latest/UCDlist.html).
        If no matching UCD type exists, this can instead be "custom:xxx",
        where xxx is an arbitrary string.  Alternatively, if the physical
        type is unknown/undefined, an element can be `None`.
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def world_axis_units(self):
        """
        Returns an iterable of strings given the units of the world
        coordinates for each axis. The strings should follow the recommended
        VOUnit standard (though as noted in the VOUnit specification
        document, units that do not follow this standard are still allowed,
        but just not recommended).
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def axis_correlation_matrix(self):
        """
        Returns an (n_world, n_pixel) matrix that indicates using booleans
        whether a given world coordinate depends on a given pixel coordinate.
        This should default to a matrix where all elements are True in the
        absence of any further information. For completely independent axes,
        the diagonal would be True and all other entries False.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def pixel_to_world_values(self, *pixel_arrays):
        """
        Convert pixel coordinates to world coordinates. This method takes
        n_pixel scalars or arrays as input, and pixel coordinates should be
        zero-based. Returns n_world scalars or arrays in units given by
        ``world_axis_units``. Note that pixel coordinates are assumed
        to be 0 at the center of the first pixel in each dimension.
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def world_to_pixel_values(self, *world_arrays):
        """
        Convert world coordinates to pixel coordinates. This method takes
        n_world scalars or arrays as input in units given by ``world_axis_units``.
        Returns n_pixel scalars or arrays. Note that pixel coordinates are
        assumed to be 0 at the center of the first pixel in each dimension.
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def world_axis_object_components(self):
        """
        A list with n_dim_world elements, where each element is a tuple with
        two items:

        * The first is a name for the world object this world array
          corresponds to, which *must* match the string names used in
          ``world_axis_object_classes``. Note that names might appear twice
          because two world arrays might correspond to a single world object
          (e.g. a celestial coordinate might have both “ra” and “dec”
          arrays, which correspond to a single sky coordinate object).

        * The second element is either a string keyword argument name or a
          positional index for the corresponding class from
          ``world_axis_object_classes``

        See below for an example of this property.
        """
        raise NotImplementedError()

    @abc.abstractproperty
    def world_axis_object_classes(self):
        """
        A dictionary with each key being a string key from
        ``world_axis_object_components``, and each value being a tuple with
        two elements:

        * The first element of the tuple must be a string specifying the
          fully-qualified name of a class, which will specify the actual
          Python object to be created.

        * The second tuple element must be a
          dictionary with the keyword arguments required to initialize the
          class.

        See below for an example of this property. Note that we don't
        require the classes to be Astropy classes since there is no
        guarantee that Astropy will have all the classes to represent all
        kinds of world coordinates. Furthermore, we recommend that the
        output be kept as human-readable as possible.
        """
        raise NotImplementedError()
