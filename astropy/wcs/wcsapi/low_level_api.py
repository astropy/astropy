import abc
import os

import numpy as np

__all__ = ["BaseLowLevelWCS", "validate_physical_types"]


class BaseLowLevelWCS(metaclass=abc.ABCMeta):
    """
    Abstract base class for the low-level WCS interface.

    This is described in `APE 14: A shared Python interface for World Coordinate
    Systems <https://doi.org/10.5281/zenodo.1188875>`_.
    """

    @property
    @abc.abstractmethod
    def pixel_n_dim(self):
        """
        The number of axes in the pixel coordinate system.
        """

    @property
    @abc.abstractmethod
    def world_n_dim(self):
        """
        The number of axes in the world coordinate system.
        """

    @property
    @abc.abstractmethod
    def world_axis_physical_types(self):
        """
        An iterable of strings describing the physical type for each world axis.

        These should be names from the VO UCD1+ controlled Vocabulary
        (http://www.ivoa.net/documents/latest/UCDlist.html). If no matching UCD
        type exists, this can instead be ``"custom:xxx"``, where ``xxx`` is an
        arbitrary string.  Alternatively, if the physical type is
        unknown/undefined, an element can be `None`.
        """

    @property
    @abc.abstractmethod
    def world_axis_units(self):
        """
        An iterable of strings given the units of the world coordinates for each
        axis.

        The strings should follow the `IVOA VOUnit standard
        <http://ivoa.net/documents/VOUnits/>`_ (though as noted in the VOUnit
        specification document, units that do not follow this standard are still
        allowed, but just not recommended).
        """

    @abc.abstractmethod
    def pixel_to_world_values(self, *pixel_arrays):
        """
        Convert pixel coordinates to world coordinates.

        This method takes `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` scalars or arrays as
        input, and pixel coordinates should be zero-based. Returns
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim` scalars or arrays in units given by
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_units`. Note that pixel coordinates are
        assumed to be 0 at the center of the first pixel in each dimension. If a
        pixel is in a region where the WCS is not defined, NaN can be returned.
        The coordinates should be specified in the ``(x, y)`` order, where for
        an image, ``x`` is the horizontal coordinate and ``y`` is the vertical
        coordinate.

        If `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim` is ``1``, this
        method returns a single scalar or array, otherwise a tuple of scalars or
        arrays is returned.
        """

    def array_index_to_world_values(self, *index_arrays):
        """
        Convert array indices to world coordinates.

        This is the same as `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values` except that
        the indices should be given in ``(i, j)`` order, where for an image
        ``i`` is the row and ``j`` is the column (i.e. the opposite order to
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values`).

        If `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim` is ``1``, this
        method returns a single scalar or array, otherwise a tuple of scalars or
        arrays is returned.
        """
        return self.pixel_to_world_values(*index_arrays[::-1])

    @abc.abstractmethod
    def world_to_pixel_values(self, *world_arrays):
        """
        Convert world coordinates to pixel coordinates.

        This method takes `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim` scalars or arrays as
        input in units given by `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_units`. Returns
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` scalars or arrays. Note that pixel
        coordinates are assumed to be 0 at the center of the first pixel in each
        dimension. If a world coordinate does not have a matching pixel
        coordinate, NaN can be returned.  The coordinates should be returned in
        the ``(x, y)`` order, where for an image, ``x`` is the horizontal
        coordinate and ``y`` is the vertical coordinate.

        If `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` is ``1``, this
        method returns a single scalar or array, otherwise a tuple of scalars or
        arrays is returned.
        """

    def world_to_array_index_values(self, *world_arrays):
        """
        Convert world coordinates to array indices.

        This is the same as `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_to_pixel_values` except that
        the indices should be returned in ``(i, j)`` order, where for an image
        ``i`` is the row and ``j`` is the column (i.e. the opposite order to
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_to_world_values`). The indices should be
        returned as rounded integers.

        If `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` is ``1``, this
        method returns a single scalar or array, otherwise a tuple of scalars or
        arrays is returned.
        """
        pixel_arrays = self.world_to_pixel_values(*world_arrays)
        if self.pixel_n_dim == 1:
            pixel_arrays = (pixel_arrays,)
        else:
            pixel_arrays = pixel_arrays[::-1]
        array_indices = tuple(
            np.asarray(np.floor(pixel + 0.5), dtype=int) for pixel in pixel_arrays
        )
        return array_indices[0] if self.pixel_n_dim == 1 else array_indices

    @property
    @abc.abstractmethod
    def world_axis_object_components(self):
        """
        A list with `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim` elements giving information
        on constructing high-level objects for the world coordinates.

        Each element of the list is a tuple with three items:

        * The first is a name for the world object this world array
          corresponds to, which *must* match the string names used in
          `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_object_classes`. Note that names might
          appear twice because two world arrays might correspond to a single
          world object (e.g. a celestial coordinate might have both “ra” and
          “dec” arrays, which correspond to a single sky coordinate object).

        * The second element is either a string keyword argument name or a
          positional index for the corresponding class from
          `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_object_classes`.

        * The third argument is a string giving the name of the property
          to access on the corresponding class from
          `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_object_classes` in
          order to get numerical values. Alternatively, this argument can be a
          callable Python object that takes a high-level coordinate object and
          returns the numerical values suitable for passing to the low-level
          WCS transformation methods.

        See the document
        `APE 14: A shared Python interface for World Coordinate Systems
        <https://doi.org/10.5281/zenodo.1188875>`_ for examples.
        """

    @property
    @abc.abstractmethod
    def world_axis_object_classes(self):
        """
        A dictionary giving information on constructing high-level objects for
        the world coordinates.

        Each key of the dictionary is a string key from
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_object_components`, and each value is a
        tuple with three elements or four elements:

        * The first element of the tuple must be a class or a string specifying
          the fully-qualified name of a class, which will specify the actual
          Python object to be created.

        * The second element, should be a tuple specifying the positional
          arguments required to initialize the class. If
          `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_object_components` specifies that the
          world coordinates should be passed as a positional argument, this this
          tuple should include `None` placeholders for the world coordinates.

        * The third tuple element must be a dictionary with the keyword
          arguments required to initialize the class.

        * Optionally, for advanced use cases, the fourth element (if present)
          should be a callable Python object that gets called instead of the
          class and gets passed the positional and keyword arguments. It should
          return an object of the type of the first element in the tuple.

        Note that we don't require the classes to be Astropy classes since there
        is no guarantee that Astropy will have all the classes to represent all
        kinds of world coordinates. Furthermore, we recommend that the output be
        kept as human-readable as possible.

        The classes used here should have the ability to do conversions by
        passing an instance as the first argument to the same class with
        different arguments (e.g. ``Time(Time(...), scale='tai')``). This is
        a requirement for the implementation of the high-level interface.

        The second and third tuple elements for each value of this dictionary
        can in turn contain either instances of classes, or if necessary can
        contain serialized versions that should take the same form as the main
        classes described above (a tuple with three elements with the fully
        qualified name of the class, then the positional arguments and the
        keyword arguments). For low-level API objects implemented in Python, we
        recommend simply returning the actual objects (not the serialized form)
        for optimal performance. Implementations should either always or never
        use serialized classes to represent Python objects, and should indicate
        which of these they follow using the
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.serialized_classes` attribute.

        See the document
        `APE 14: A shared Python interface for World Coordinate Systems
        <https://doi.org/10.5281/zenodo.1188875>`_ for examples .
        """

    # The following three properties have default fallback implementations, so
    # they are not abstract.

    @property
    def array_shape(self):
        """
        The shape of the data that the WCS applies to as a tuple of length
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` in ``(row, column)``
        order (the convention for arrays in Python).

        If the WCS is valid in the context of a dataset with a particular
        shape, then this property can be used to store the shape of the
        data. This can be used for example if implementing slicing of WCS
        objects. This is an optional property, and it should return `None`
        if a shape is not known or relevant.
        """
        if self.pixel_shape is None:
            return None
        else:
            return self.pixel_shape[::-1]

    @property
    def pixel_shape(self):
        """
        The shape of the data that the WCS applies to as a tuple of length
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim` in ``(x, y)``
        order (where for an image, ``x`` is the horizontal coordinate and ``y``
        is the vertical coordinate).

        If the WCS is valid in the context of a dataset with a particular
        shape, then this property can be used to store the shape of the
        data. This can be used for example if implementing slicing of WCS
        objects. This is an optional property, and it should return `None`
        if a shape is not known or relevant.

        If you are interested in getting a shape that is comparable to that of
        a Numpy array, you should use
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.array_shape` instead.
        """
        return None

    @property
    def pixel_bounds(self):
        """
        The bounds (in pixel coordinates) inside which the WCS is defined,
        as a list with `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim`
        ``(min, max)`` tuples.

        The bounds should be given in ``[(xmin, xmax), (ymin, ymax)]``
        order. WCS solutions are sometimes only guaranteed to be accurate
        within a certain range of pixel values, for example when defining a
        WCS that includes fitted distortions. This is an optional property,
        and it should return `None` if a shape is not known or relevant.
        """
        return None

    @property
    def pixel_axis_names(self):
        """
        An iterable of strings describing the name for each pixel axis.

        If an axis does not have a name, an empty string should be returned
        (this is the default behavior for all axes if a subclass does not
        override this property). Note that these names are just for display
        purposes and are not standardized.
        """
        return [""] * self.pixel_n_dim

    @property
    def world_axis_names(self):
        """
        An iterable of strings describing the name for each world axis.

        If an axis does not have a name, an empty string should be returned
        (this is the default behavior for all axes if a subclass does not
        override this property). Note that these names are just for display
        purposes and are not standardized. For standardized axis types, see
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.world_axis_physical_types`.
        """
        return [""] * self.world_n_dim

    @property
    def axis_correlation_matrix(self):
        """
        Returns an (`~astropy.wcs.wcsapi.BaseLowLevelWCS.world_n_dim`,
        `~astropy.wcs.wcsapi.BaseLowLevelWCS.pixel_n_dim`) matrix that
        indicates using booleans whether a given world coordinate depends on a
        given pixel coordinate.

        This defaults to a matrix where all elements are `True` in the absence
        of any further information. For completely independent axes, the
        diagonal would be `True` and all other entries `False`.
        """
        return np.ones((self.world_n_dim, self.pixel_n_dim), dtype=bool)

    @property
    def serialized_classes(self):
        """
        Indicates whether Python objects are given in serialized form or as
        actual Python objects.
        """
        return False

    def _as_mpl_axes(self):
        """Compatibility hook for Matplotlib and WCSAxes.

        With this method, one can do::

            from astropy.wcs import WCS
            import matplotlib.pyplot as plt
            wcs = WCS('filename.fits')
            fig = plt.figure()
            ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
            ...

        and this will generate a plot with the correct WCS coordinates on the
        axes.
        """
        from astropy.visualization.wcsaxes import WCSAxes

        return WCSAxes, {"wcs": self}


UCDS_FILE = os.path.join(os.path.dirname(__file__), "data", "ucds.txt")
with open(UCDS_FILE) as f:
    VALID_UCDS = {x.strip() for x in f.read().splitlines()[1:]}


def validate_physical_types(physical_types):
    """
    Validate a list of physical types against the UCD1+ standard.
    """
    for physical_type in physical_types:
        if (
            physical_type is not None
            and physical_type not in VALID_UCDS
            and not physical_type.startswith("custom:")
        ):
            raise ValueError(
                f"'{physical_type}' is not a valid IOVA UCD1+ physical type. It must be"
                " a string specified in the list"
                " (http://www.ivoa.net/documents/latest/UCDlist.html) or if no"
                " matching type exists it can be any string prepended with 'custom:'."
            )
