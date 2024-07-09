# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Defines a class used to store and interface with covariance matrices.
"""

import warnings
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from scipy import sparse

from astropy import log
from astropy.io import fits

__all__ = ["Covariance"]


class Covariance:
    r"""
    A general utility for storing, manipulating, and file I/O of sparse
    covariance matrices.

    .. warning::

        Although possible to abstract further, use of ``raw_shape``
        can only be 2D in the current implementation.

    Parameters
    ----------
    inp : `scipy.sparse.csr_matrix`, optional
        Covariance matrix to store. Input **must** be covariance data, not
        correlation data.  If None, the covariance object is instantiated as
        being empty.

    impose_triu : :obj:`bool`, optional
        Flag to force the `scipy.sparse.csr_matrix` object to only be the upper
        triangle of the covariance matrix. Covariance matrices are symmetric by
        definition, :math:`C_{ij} = C_{ji}`, so it's not necessary to keep both
        values. This flag will force a call to `scipy.sparse.triu` when setting
        the covariance matrix.  If False, the input matrix is
        **assumed** to only have the upper triangle of numbers.

    correlation : :obj:`bool`, optional
        Convert the input to a correlation matrix. The input
        **must** be the covariance matrix.

    raw_shape : :obj:`tuple`, optional
        The covariance data is for a higher dimensional array with this shape.
        For example, if the covariance data is for a 2D image with shape
        ``(nx,ny)`` -- the shape of the covariance array is be ``(nx*ny,
        nx*ny)`` -- set ``raw_shape=(nx,ny)``. This is primarily used for
        reading and writing; see also :func:`transpose_raw_shape`.  If None, any
        higher dimensionality is ignored.  Currently, this is only allowed if
        the raw shape is 2D.

    Raises
    ------
    TypeError
        Raised if the input array not a `scipy.sparse.csr_matrix` object.

    NotImplementedError
        Raised if ``raw_shape`` is not a two-tuple.

    ValueError
        Raised if ``raw_shape`` is provided and the input covariance matrix
        ``inp`` does not have the expected shape.  I.e., if
        ``raw_shape=(nx,ny)``, the covariance matrix must have the shape
        ``(nx*ny,nx*ny)``.

    Attributes
    ----------
    cov : `scipy.sparse.csr_matrix`
        The covariance matrix stored in sparse format.

    shape : :obj:`tuple`
        Shape of the full array.

    raw_shape : :obj:`tuple`
        The covariance data is for a higher dimensional array with this shape.
        For example, if the covariance data is for a 2D image, this would be
        ``(nx,ny)`` and the shape of the covariance array would be ``(nx*ny,
        nx*ny)``.

    nnz : :obj:`int`
        The number of non-zero covariance matrix elements.

    var : `numpy.ndarray`
        Array with the variance provided by the diagonal of the covariance
        matrix. This is only populated if necessary, either by being requested
        (:func:`variance`) or if needed to convert between covariance and
        correlation matrices.

    is_correlation : :obj:`bool`
        Flag that the covariance matrix has been saved as a variance vector and
        a correlation matrix.
    """

    def __init__(self, inp, impose_triu=False, correlation=False, raw_shape=None):
        self.cov = inp
        self.shape = None
        self.raw_shape = raw_shape
        if self.raw_shape is not None and self.raw_shape.ndim != 2:
            raise NotImplementedError("The source raw_shape can only be 2D.")

        self.nnz = None
        self.var = None
        self.is_correlation = False

        # Return empty object
        if self.cov is None:
            return

        # Set the dimensionality, check that each element of the array
        # has the correct type, count the number of non-zero covariance
        # values, and (if necessary) initialize the indices for each
        # covariance matrix
        if not sparse.isspmatrix_csr(self.cov):
            raise TypeError(
                "Input covariance matrix (or elements) must be csr matrices."
            )
        self.nnz = self.cov.nnz

        # Set the shape of the full matrix/cube
        self.shape = self.cov.shape
        if self.raw_shape is not None and np.prod(self.raw_shape) != self.shape[0]:
            raise ValueError(
                "Product of raw shape must match the covariance axis length."
            )

        # If requested, impose that the input matrix only have values in
        # its upper triangle.
        if impose_triu:
            self._impose_upper_triangle()

        # Set the variance array and the correlation matrix flag
        if correlation:
            self.to_correlation()

    @classmethod
    def from_samples(cls, samples, cov_tol=None, rho_tol=None):
        r"""
        Build a covariance object using discrete samples.

        The covariance is generated using `numpy.cov` for a set of discretely
        sampled data for an :math:`N`-dimensional parameter space.

        Parameters
        ----------
        samples : `numpy.ndarray`
            Array with samples drawn from an :math:`N`-dimensional parameter
            space. The shape of the input array must be :math:`N_{\rm par}\times
            N_{\rm samples}`.

        cov_tol : :obj:`float`, optional
            Any covariance value less than this is assumed to be equivalent to
            (and set to) 0.

        rho_tol : :obj:`float`, optional
            Any correlation coefficient less than this is assumed to be
            equivalent to (and set to) 0.

        Returns
        -------
        :class:`Covariance`
            An :math:`N_{\rm par}\times N_{\rm par}` covariance matrix built
            using the provided samples.

        Raises
        ------
        ValueError
            Raised if the input array is not 2D or if the number of samples (length
            of the second axis) is less than 2.
        """
        if samples.ndim != 2:
            raise ValueError("Input samples for covariance matrix must be a 2D array!")
        if samples.shape[1] < 2:
            raise ValueError("Fewer than two samples provided!")
        return Covariance.from_array(np.cov(samples), cov_tol=cov_tol, rho_tol=rho_tol)

    @classmethod
    def from_array(cls, covar, cov_tol=None, rho_tol=None, raw_shape=None):
        r"""
        Define a covariance object using a dense array.

        Parameters
        ----------
        covar : array-like
            Array with the covariance data. The shape of the array must be
            square. Input can be any object that can be converted to a dense
            array using the object method ``toarray`` or using
            ``numpy.atleast_2d``.

        cov_tol : :obj:`float`, optional
            Any covariance value less than this is assumed to be equivalent to
            (and set to) 0.

        rho_tol : :obj:`float`, optional
            Any correlation coefficient less than this is assumed to be
            equivalent to (and set to) 0.

        raw_shape : :obj:`tuple`, optional
            The covariance data is for a higher dimensional array with this
            shape.  For example, if the covariance data is for a 2D image with
            shape ``(nx,ny)`` -- the shape of the covariance array is be
            ``(nx*ny, nx*ny)`` -- set ``raw_shape=(nx,ny)``. This is primarily
            used for reading and writing; see also :func:`transpose_raw_shape`.
            If None, any higher dimensionality is ignored.  Currently, this is
            only allowed if the raw shape is 2D.

        Returns
        -------
        :class:`Covariance`
            The covariance matrix built using the provided array.

        Raises
        ------
        ValueError
            Raised if ``covar`` could not be converted to a dense array.
        """
        try:
            _covar = covar.toarray()
        except AttributeError:
            _covar = np.atleast_2d(covar)
        if not isinstance(_covar, np.ndarray) or _covar.ndim != 2:
            raise ValueError(
                "Could not convert input covariance data into a 2D dense array."
            )

        n = _covar.shape[0]
        if rho_tol is not None:
            variance = np.diag(_covar)
            correlation = _covar / np.ma.sqrt(variance[:, None] * variance[None, :])
            correlation[np.ma.absolute(correlation) < rho_tol] = 0.0
            _covar = correlation.filled(0.0) * np.ma.sqrt(
                variance[:, None] * variance[None, :]
            ).filled(0.0)
        if cov_tol is not None:
            _covar[_covar < cov_tol] = 0.0

        indx = _covar > 0.0
        i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        return cls(
            sparse.coo_matrix(
                (_covar[indx].ravel(), (i[indx].ravel(), j[indx].ravel())), shape=(n, n)
            ).tocsr(),
            impose_triu=True,
            raw_shape=raw_shape,
        )

    @classmethod
    def from_fits(
        cls,
        source,
        ivar_ext="IVAR",
        transpose_ivar=False,
        covar_ext="CORREL",
        correlation=False,
        quiet=False,
    ):
        r"""
        Read covariance data from a FITS file.

        This read operation matches the data saved to a FITS file using
        :func:`write`. The class can read covariance data written by other
        programs *as long as they have a commensurate format*. See the
        description of the :func:`write` method.

        If the extension names and column names are correct, :class:`Covariance`
        can read FITS files that were not produced explicitly by this method.
        This is useful for files that include the covariance data extensions
        among others. The methods :func:`output_hdus` and
        :func:`coordinate_data` are provided to produce the data that can be
        placed in any FITS file.

        The method determines if the output data were reshaped by checking the
        number of columns in the binary table.

        When the covariance data are for a higher dimensional array, the memory
        order of the higher dimensional array (particularly whether it's
        constructed row- or column-major) is important to ensuring the loaded
        data is correct. This is why we have chosen the specific format for the
        binary table used to store the covariance data. Regardless of whether
        the data was written by a row-major or column-major language, the format
        should be such that this class can properly read and recover the
        covariance matrix. However, in some cases, you still may need to
        transpose the inverse variance data; see ``transpose_ivar``.

        Parameters
        ----------
        source : :obj:`str`, `Path`, `astropy.io.fits.HDUList`
            Initialize the object using an `astropy.io.fits.HDUList` object or
            path to a FITS file.

        ivar_ext : :obj:`int`, :obj:`str`, optional
            The name or index of the extension (see ``source``) with the inverse
            variance data.  If None, the variance is taken as unity.

        transpose_ivar : :obj:`bool`, optional
            Flag to transpose the inverse variance data before rescaling the
            correlation matrix. Should only be necessary in some cases when the
            covariance data was written by a method other than :func:`write`.

        covar_ext : :obj:`int`, :obj:`str`, optional
            The name or index of the extension with covariance data.  This
            **cannot** be None.

        correlation : :obj:`bool`, optional
            Return the matrix as a correlation matrix.  If False, use the data
            (always saved in correlation format; see :func:`write`) to construct
            the covariance matrix.

        quiet : :obj:`bool`, optional
            Suppress terminal output.

        Returns
        -------
        :class:`Covariance`
            The covariance matrix read from the provided source.

        Raises
        ------
        ValueError
            Raised if ``covar_ext`` is None.
        """
        # Check input
        if covar_ext is None:
            raise ValueError(
                "Must provide extension name or index with the covariance data."
            )

        # Open the provided source, if it hasn't been yet
        source_is_hdu = isinstance(source, fits.HDUList)
        hdu = source if source_is_hdu else fits.open(source)

        # Read coordinate data
        shape = eval(hdu[covar_ext].header["COVSHAPE"])
        raw_shape = (
            eval(hdu[covar_ext].header["COVRWSHP"])
            if "COVRWSHP" in hdu[covar_ext].header
            else None
        )
        reshape = len(hdu[covar_ext].columns.names) in [5, 6]
        if reshape:
            i_c1, i_c2, j_c1, j_c2, rhoij = (
                hdu[covar_ext].data[ext].copy()
                for ext in ["INDXI_C1", "INDXI_C2", "INDXJ_C1", "INDXJ_C2", "RHOIJ"]
            )
            if raw_shape is None:
                raw_shape = Covariance.square_shape(shape[0])
            i = np.ravel_multi_index((i_c1, i_c2), raw_shape)
            j = np.ravel_multi_index((j_c1, j_c2), raw_shape)
            # Make sure the data are only in the upper triangle
            indx = j < i
            i[indx], j[indx] = j[indx], i[indx]
        else:
            i, j, rhoij = (
                hdu[covar_ext].data[ext] for ext in ["INDXI", "INDXJ", "RHOIJ"]
            )

        # Number of non-zero elements
        nnz = len(rhoij)

        # Inverse variance data
        ivar = None if ivar_ext is None else hdu[ivar_ext].data

        # Done with the hdu so close it, if necessary
        if not source_is_hdu:
            hdu.close()

        # Set correlation data
        if ivar is not None:
            if transpose_ivar:
                ivar = ivar.T
            if reshape:
                ivar = ivar.ravel()
        var = (
            np.ones(shape[1:], dtype=float)
            if ivar is None
            else np.ma.power(ivar, -1).filled(0.0)
        )
        cij = rhoij * np.sqrt(var[i] * var[j])
        cov = sparse.coo_matrix((cij, (i, j)), shape=shape).tocsr()

        # Report
        if not quiet:
            log.info("Read covariance cube:")
            log.info(
                f'       output type: {"Correlation" if correlation else "Covariance"}'
            )
            log.info(f"             shape: {shape}")
            log.info(f"   non-zero values: {nnz}")

        return cls(cov, correlation=correlation, raw_shape=raw_shape)

    @classmethod
    def from_matrix_multiplication(cls, T, Sigma):
        r"""
        Construct the covariance matrix that results from a matrix
        multiplication.

        The matrix multiplication should be of the form:

        .. math::

            {\mathbf T} \times {\mathbf X} = {\mathbf Y}

        where :math:`{\mathbf T}` is a transfer matrix of size :math:`N_y\times
        N_x`, :math:`{\mathbf X}` is a vector of size :math:`N_x`, and
        :math:`{\mathbf Y}` is the vector of length :math:`{N_y}` that results
        from the multiplication.

        The covariance matrix is then

        .. math::

             {\mathbf C} = {\mathbf T} \times {\mathbf \Sigma} \times
             {\mathbf T}^{\rm T},

        where :math:`{\mathbf \Sigma}` is the covariance matrix for the elements
        of :math:`{\mathbf X}`. If ``Sigma`` is provided as a vector of length
        :math:`N_x`, it is assumed that the elements of :math:`{\mathbf X}` are
        independent and the provided vector gives the *variance* in each
        element; i.e., the provided data represent the diagonal of
        :math:`{\mathbf \Sigma}`.

        Parameters
        ----------
        T : `scipy.sparse.csr_matrix`, `numpy.ndarray`
            Transfer matrix.  See above.

        Sigma : `scipy.sparse.csr_matrix`, `numpy.ndarray`
            Covariance matrix.  See above.

        Returns
        -------
        :class:`Covariance`
            The covariance matrix resulting from the matrix multiplication.

        Raises
        ------
        ValueError
            Raised if the provided arrays are not two dimensional or if there is
            a shape mismatch.
        """
        if T.ndim != 2:
            raise ValueError("Input transfer matrix must be two-dimensional.")
        nx = T.shape[1]
        if Sigma.shape != (nx, nx) and Sigma.shape != (nx,):
            raise ValueError(
                "Shape of input variance matrix must be either "
                f"({nx},{nx}) or ({nx},)."
            )
        # If it isn't already, convert T to a csr_matrix
        _T = T if isinstance(T, sparse.csr_matrix) else sparse.csr_matrix(T)
        # Set the covariance matrix in X
        _Sigma = (
            sparse.coo_matrix(
                (Sigma, (np.arange(nx), np.arange(nx))), shape=(nx, nx)
            ).tocsr()
            if Sigma.ndim == 1
            else (
                Sigma
                if isinstance(Sigma, sparse.csr_matrix)
                else sparse.csr_matrix(Sigma)
            )
        )
        # Construct the covariance matrix
        return cls(sparse.triu(_T.dot(_Sigma.dot(_T.transpose()))).tocsr())

    @classmethod
    def from_variance(cls, variance, correlation=False):
        r"""
        Construct a diagonal covariance matrix using the provided variance.

        Parameters
        ----------
        variance : `numpy.ndarray`
            The variance vector.

        correlation : :obj:`bool`, optional
            Upon instantiation, convert the :class:`Covariance` object to a
            correlation matrix.

        Returns
        -------
        :class:`Covariance`
            The diagonal covariance matrix.
        """
        return cls(sparse.csr_matrix(np.diagflat(variance)), correlation=correlation)

    def _impose_upper_triangle(self):
        """
        Force :attr:`cov` to only contain non-zero elements in its upper
        triangle.
        """
        self.cov = sparse.triu(self.cov).tocsr()
        self.nnz = self.cov.nnz

    def full(self):
        r"""
        Return a `scipy.sparse.csr_matrix` object with both its upper and lower
        triangle filled, ensuring that they are symmetric.

        This method is essentially equivalent to :func:`toarray`
        except that it returns a sparse array.

        Returns
        -------
        `scipy.sparse.csr_matrix`
            The sparse matrix with both the upper and lower triangles filled
            (with symmetric information).
        """
        a = self.cov
        return sparse.triu(a) + sparse.triu(a, 1).T

    @staticmethod
    def square_shape(size):
        """
        Determine the shape of a 2D square array resulting from reshaping a
        vector.

        Parameters
        ----------
        size : :obj:`int`
            Size of the vector.

        Returns
        -------
        :obj:`int`
            Length of one axis of the square array.

        Raises
        ------
        ValueError
            Raised if the provided ``size`` is not the square of an integer.
        """
        # Get the size of the square image (and make sure it's square)
        n = np.floor(np.sqrt(size)).astype(int)
        if n * n != size:
            raise ValueError(f"{size} is not the square of an integer!")
        return (n, n)

    def transpose_raw_shape(self):
        """
        Reorder the covariance array to account for a transpose in the ravel
        ordering of the source array.

        For a higher dimensional source array (see :attr:`raw_shape`), the
        indices in the covariance matrix relevant to source array can be found
        using `numpy.ravel_multi_index`_ (or `numpy.unravel_index`_ for vice
        versa). However, if the source array is transposed, these operations
        will be invalid.

        This method constructs a new :class:`Covariance` object from the
        existing data but with the indices rearranged for the transposed source
        array.

        .. warning::

            If :attr:`raw_shape` it is not defined, a warning is
            issued and this method simply returns a copy of the
            current covariance data.

        Returns
        -------
        :class:`Covariance`
            The covariance matrix with a transposed ordering of the raw shape.
        """
        if self.raw_shape is None:
            warnings.warn("Covariance array raw shape undefined.  Returning a copy.")
            return self.copy()

        raw_shape_t = self.raw_shape[::-1]

        # Get the reshaped coordinate data
        i_c1, i_c2, j_c1, j_c2, rhoij, var = self.coordinate_data(reshape=True)

        # Get the current covariance data. TODO: Allow coordinate
        # data to return covariance data instead of it always being
        # correlation data.
        i = np.ravel_multi_index((i_c1, i_c2), self.raw_shape)
        j = np.ravel_multi_index((j_c1, j_c2), self.raw_shape)
        cij = rhoij * np.sqrt(var.flat[i] * var.flat[j])

        # Get the new coordinates
        i = np.ravel_multi_index((i_c2, i_c1), raw_shape_t)
        j = np.ravel_multi_index((j_c2, j_c1), raw_shape_t)
        indx = j < i
        i[indx], j[indx] = j[indx], i[indx]

        # Return the new covariance matrix
        return Covariance(
            sparse.coo_matrix((cij, (i, j)), shape=self.shape).tocsr(),
            raw_shape=raw_shape_t,
        )

    def apply_new_variance(self, var):
        """
        Using the same correlation coefficients, return a new
        :class:`Covariance` object with the provided variance.

        Parameters
        ----------
        var : `numpy.ndarray`
            Variance vector. Must have a length that matches the shape of this
            :class:`Covariance` instance.

        Returns
        -------
        :class:`Covariance`
            A covariance matrix with the same shape and correlation coefficients
            and this object, but with the provided variance.

        Raises
        ------
        ValueError
            Raised if the length of the variance vector is incorrect.
        """
        if var.shape != self.shape[1:]:
            raise ValueError(
                f"Provided variance has incorrect shape.  Expected {self.shape[1:]}, "
                f"provided {var.shape}."
            )

        # Convert to a correlation matrix, if needed
        is_correlation = self.is_correlation
        if not is_correlation:
            self.to_correlation()

        # Pull out the non-zero values
        i, j, c = sparse.find(self.cov)
        # Apply the new variance
        new_cov = sparse.coo_matrix(
            (c * np.sqrt(var[i] * var[j]), (i, j)), shape=self.shape
        ).tocsr()

        # Revert to covariance matrix, if needed
        if not is_correlation:
            self.revert_correlation()

        # Return a new covariance matrix
        return Covariance(new_cov, correlation=is_correlation)

    def copy(self):
        """
        Return a copy of this Covariance object.

        Returns
        -------
        :class:`Covariance`
            A copy of the current covariance matrix.
        """
        # If the data is saved as a correlation matrix, first revert to
        # a covariance matrix
        is_correlation = self.is_correlation
        if self.is_correlation:
            self.revert_correlation()

        # Create the new Covariance instance with a copy of the data
        cp = Covariance(self.cov.copy())

        # If necessary, convert the data to a correlation matrix
        if is_correlation:
            self.to_correlation()
            cp.to_correlation()
        return cp

    def toarray(self):
        """
        Convert the sparse covariance matrix to a dense array, filled
        with zeros where appropriate.

        Returns
        -------
        `numpy.ndarray`
            Dense array with the full covariance matrix.
        """
        return self.full().toarray()

    def show(self, zoom=None, ofile=None, log10=False):
        """
        Show a covariance/correlation matrix data.

        This converts the covariance matrix to a filled array and plots the
        array using `matplotlib.pyplot.imshow`. If an output file is provided,
        the image is redirected to the designated output file; otherwise, the
        image is plotted to the screen.

        Parameters
        ----------
        zoom : :obj:`float`, optional
            Factor by which to zoom in on the center of the image by *removing
            the other regions of the array*. E.g. ``zoom=2`` will show only the
            central quarter of the covariance matrix.

        ofile : :obj:`str`, optional
            If provided, the plot is output to this file instead of being
            plotted to the screen.

        log10 : :obj:`bool`, optional
            Plot the base-10 log of the covariance value.
        """
        # Convert the covariance matrix to an array
        a = self.toarray()

        # Remove some fraction of the array to effectively zoom in on
        # the center of the covariance matrix
        if zoom is not None:
            xs = int(self.shape[0] / 2 - self.shape[0] / 2 / zoom)
            xe = xs + int(self.shape[0] / zoom) + 1
            ys = int(self.shape[1] / 2 - self.shape[1] / 2 / zoom)
            ye = ys + int(self.shape[1] / zoom) + 1
            a = a[xs:xe, ys:ye]

        # Create the plot
        fig = plt.figure(1)
        im = plt.imshow(
            np.ma.log10(a) if log10 else a, interpolation="nearest", origin="lower"
        )
        plt.colorbar()
        if ofile is None:
            # Print the plot to the screen if no output file is provided.
            plt.show()
        else:
            # Print the plot to the designated output file
            fig.canvas.print_figure(ofile)
        fig.clear()
        plt.close(fig)

    # TODO: Should there be any difference between this and `coordinate_data`
    def find(self):
        """
        Find the non-zero values in the **full** covariance matrix (not just the
        upper triangle).

        This is a simple wrapper for :func:`full` and `scipy.sparse.find`.

        Returns
        -------
        tuple
            A tuple of arrays ``i``, ``j``, and ``c``. The arrays ``i`` and
            ``j`` contain the index coordinates of the non-zero values, and
            ``c`` contains the values themselves.
        """
        return sparse.find(self.full())

    def coordinate_data(self, reshape=False):
        r"""
        Construct data arrays with the non-zero covariance components in
        coordinate format.

        This procedure is primarily used when constructing the data arrays to
        write to a FITS file.

        Regardless of whether or not the current internal data is a covariance
        matrix or a correlation matrix, the data is always returned as a
        correlation matrix.

        Matching the class convention, the returned data only includes the upper
        triangle.

        Four columns are returned that contain:

            - :math:`i,j`: The row and column indices, respectively, of the
              covariance matrix.

            - :math:`rho_{ij}`: The correlation coefficient between pixels
              :math:`i` and :math:`j`.

            - :math:`V_i`: The variance in each pixel; i.e., the value of
              :math:`C_{ii} \forall i`. If reshape is True, this will be output
              as a two-dimensional array.

        If using the reshape option, the :math:`i,j` indices are
        converted to two columns each providing the indices in the
        associated reshaped array with coordinates :math:`c_1,c_2`.

        Parameters
        ----------
        reshape : :obj:`bool`, optional
            Reshape the output in the case when :math:`i,j` are actually pixels
            in a two-dimensional image. The shape of the image is expected to be
            square, such that the shape of the covariance matrix is
            :math:`N_x\times N_y`, where :math:`N_x = N_y`. Each of the ``I``
            and ``J`` output columns are then split into two columns according
            to associate coordinates, such that there are 8 output columns.

        Returns
        -------
        :obj:`tuple`

            Four or six `numpy.ndarray` objects depending on if the data has
            been reshaped into an image.  In detail:

                - If not reshaping: The four returned objects contain the
                  indices along the first and second axes, the correlation
                  coefficient, and the variance vector.

                - If reshaping: The six returned objects are the unraveled
                  indices in the 2D source array (two indices each for the two
                  axes of the covariance matrix, ordered as the first and second
                  indices for the first covariance index then for the second
                  covariance index), the correlation coefficient, and the 2D
                  variance array.
        """
        # Only ever print correlation matrices
        is_correlation = self.is_correlation
        if not is_correlation:
            self.to_correlation()

        if reshape:
            # If reshaping, get the new shape; assume the 2D array is
            # square if raw_shape is None.
            new_shape = (
                Covariance.square_shape(self.shape[0])
                if self.raw_shape is None
                else self.raw_shape
            )

        # Get the data
        i, j, rhoij = sparse.find(self.cov)

        # If object was originally a covariance matrix, revert it back
        if not is_correlation:
            self.revert_correlation()

        # Return the data
        if not reshape:
            # Returns four arrays
            return i, j, rhoij, self.var.copy()

        # Returns six arrays
        i_c1, i_c2 = np.unravel_index(i, new_shape)
        j_c1, j_c2 = np.unravel_index(j, new_shape)
        return i_c1, i_c2, j_c1, j_c2, rhoij, self.var.reshape(new_shape).copy()

    def output_hdus(self, reshape=False, hdr=None):
        r"""
        Construct the output HDUs and header that contain the covariance data.

        Parameters
        ----------
        reshape : :obj:`bool`, optional
            Reshape the output in the case when :math:`i,j` are actually pixels
            in a two-dimensional image. The shape of the image is expected to be
            square, such that the shape of the covariance matrix is
            :math:`N_x\times N_y`, where :math:`N_x = N_y`. Each of the ``I``
            and ``J`` output columns are then split into two columns according
            to associate coordinates, such that there are 8 output columns.

        hdr : `astropy.io.fits.Header`, optional
            `astropy.io.fits.Header` instance to which to add covariance
            keywords. If None, a new `astropy.io.fits.Header` instance is
            returned.

        Returns
        -------
        tuple
            Returns three objects: (1) The header for the primary HDU; (2) an
            `astropy.io.fits.ImageHDU` object with the variance vector; and (3)
            an `astropy.io.fits.BinTableHDU` object with the correlation
            coefficients.
        """
        # Use input header or create a minimal one
        _hdr = fits.Header() if hdr is None else hdr
        # Ensure the input header has the correct type
        if not isinstance(_hdr, fits.Header):
            raise TypeError("Input header must have type astropy.io.fits.Header.")
        _hdr["COVSHAPE"] = (str(self.shape), "Shape of the correlation matrix")

        tbl_hdr = fits.Header()
        tbl_hdr["COVSHAPE"] = (str(self.shape), "Shape of the correlation matrix")
        if self.raw_shape is not None:
            tbl_hdr["COVRWSHP"] = (str(self.raw_shape), "Raw shape of the source data")

        # Construct the correlation binary table HDU
        if reshape:
            i_c1, i_c2, j_c1, j_c2, rhoij, var = self.coordinate_data(reshape=reshape)
            covar_hdu = fits.BinTableHDU.from_columns(
                [
                    fits.Column(name="INDXI_C1", format="1J", array=i_c1),
                    fits.Column(name="INDXI_C2", format="1J", array=i_c2),
                    fits.Column(name="INDXJ_C1", format="1J", array=j_c1),
                    fits.Column(name="INDXJ_C2", format="1J", array=j_c2),
                    fits.Column(name="RHOIJ", format="1D", array=rhoij),
                ],
                name="CORREL",
                header=tbl_hdr,
            )
        else:
            i, j, rhoij, var = self.coordinate_data(reshape=reshape)
            covar_hdu = fits.BinTableHDU.from_columns(
                [
                    fits.Column(name="INDXI", format="1J", array=i),
                    fits.Column(name="INDXJ", format="1J", array=j),
                    fits.Column(name="RHOIJ", format="1D", array=rhoij),
                ],
                name="CORREL",
                header=tbl_hdr,
            )

        ivar_hdu = fits.ImageHDU(data=np.ma.power(var, -1.0).filled(0.0), name="IVAR")
        return _hdr, ivar_hdu, covar_hdu

    def write(self, ofile, reshape=False, hdr=None, overwrite=False):
        r"""
        Write the covariance object to a FITS file.

        Objects written using this function can be reinstantiated using
        :func:`from_fits`.

        The covariance matrix is stored in "coordinate" format using FITS binary
        tables; see `scipy.sparse.coo_matrix`. The matrix is *always* stored as
        a correlation matrix, even if the object is currently in the state
        holding the covariance data.

        Independent of the dimensionality of the covariance matrix, the written
        file has a ``PRIMARY`` extension with the keyword ``COVSHAPE`` that
        specifies the original dimensions of the covariance matrix; see
        :attr:`shape`.

        The correlation data are written to the ``CORREL`` extension.  The
        number of columns in this extension depends on the provided keywords;
        see :func:`coordinate_data`. The column names are:

            - ``INDXI``, ``INDXJ``: indices in the covariance matrix.  ``INDXI``
              and ``INDXJ`` are separated into two columns if the output is
              reshaped; these columns are ``INDXI_C1``, ``INDXI_C2``,
              ``INDXJ_C1``, ``INDXJ_C2``.

            - ``RHOIJ``: The non-zero correlation coefficients located
              the specified coordinates

        The inverse of the variance along the diagonal of the covariance matrix
        is output in an ImageHDU in extension ``IVAR``.

        Parameters
        ----------
        ofile : :obj:`str`
            File name for the output.

        reshape : :obj:`bool`, optional
            Reshape the output in the case when :math:`i,j` are actually pixels
            in a two-dimensional image. The shape of the image is expected to be
            square, such that the shape of the covariance matrix is
            :math:`N_x\times N_y`, where :math:`N_x = N_y`. Each of the ``I``
            and ``J`` output columns are then split into two columns according
            to associate coordinates.

        hdr : `astropy.io.fits.Header`_, optional
            A header object to include in the PRIMARY extension.  The SHAPE
            keyword will be added/overwritten.

        overwrite : :obj:`bool`, optional
            Overwrite any existing file.

        Raises
        ------
        FileExistsError
            Raised if the output file already exists and overwrite is False.

        TypeError
            Raised if the input ``hdr`` does not have the correct type.
        """
        _ofile = Path(ofile).absolute()
        if _ofile.is_file() and not overwrite:
            raise FileExistsError(
                f"{ofile} exists!  Use 'overwrite=True' to overwrite."
            )

        # Construct HDUList and write the FITS file
        _hdr, ivar_hdu, covar_hdu = self.output_hdus(reshape=reshape, hdr=hdr)
        fits.HDUList([fits.PrimaryHDU(header=_hdr), ivar_hdu, covar_hdu]).writeto(
            ofile, overwrite=overwrite, checksum=True
        )

    def variance(self, copy=True):
        """
        Return the variance vector of the covariance matrix.

        Parameters
        ----------
        copy : :obj:`bool`, optional
            Return a copy instead of a reference.
        """
        if self.var is not None:
            return self.var.copy() if copy else self.var

        self.var = np.diag(self.cov.toarray()).copy()
        return self.var

    def to_correlation(self):
        r"""
        Convert the covariance matrix into a correlation matrix by dividing each
        element by the variances.

        If the matrix is a correlation matrix already (see
        :attr:`is_correlation`), no operations are performed.  Otherwise, the
        variance vectors are computed, if necessary, and used to normalize the
        covariance values.

        A :class:`Covariance` object can be reverted from a correlation matrix
        using :func:`revert_correlation`.
        """
        # Object is already a correlation matrix
        if self.is_correlation:
            return

        # Ensure that the variance has been calculated
        self.variance()

        self.is_correlation = True
        i, j, c = sparse.find(self.cov)
        self.cov = sparse.coo_matrix(
            (c / np.sqrt(self.var[i] * self.var[j]), (i, j)), shape=self.shape
        ).tocsr()

    def revert_correlation(self):
        r"""
        Revert the object from a correlation matrix back to a full covariance
        matrix.

        This function does nothing if the correlation flag has not been flipped.
        The variances must have already been calculated!
        """
        if not self.is_correlation:
            return

        i, j, c = sparse.find(self.cov)
        self.cov = sparse.coo_matrix(
            (c * np.sqrt(self.var[i] * self.var[j]), (i, j)), shape=self.shape
        ).tocsr()
        self.is_correlation = False
