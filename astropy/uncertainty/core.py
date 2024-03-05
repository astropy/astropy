# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Distribution class and associated machinery.
"""
import builtins

import numpy as np

from astropy import stats
from astropy import units as u
from astropy.utils.compat.numpycompat import NUMPY_LT_1_23, NUMPY_LT_2_0

if NUMPY_LT_2_0:
    from numpy.core.multiarray import normalize_axis_index
    from numpy.lib.function_base import _parse_gufunc_signature
    from numpy.lib.stride_tricks import DummyArray
else:
    from numpy.lib._function_base_impl import _parse_gufunc_signature
    from numpy.lib._stride_tricks_impl import DummyArray
    from numpy.lib.array_utils import normalize_axis_index


from .function_helpers import FUNCTION_HELPERS

__all__ = ["Distribution"]


# we set this by hand because the symbolic expression (below) requires scipy
# SMAD_SCALE_FACTOR = 1 / scipy.stats.norm.ppf(0.75)
SMAD_SCALE_FACTOR = 1.48260221850560203193936104071326553821563720703125


class Distribution:
    """A scalar value or array values with associated uncertainty distribution.

    This object will take its exact type from whatever the ``samples``
    argument is. In general this is expected to be ``NdarrayDistribution`` for
    |ndarray| input, and, e.g., ``QuantityDistribution`` for a subclass such
    as |Quantity|. But anything compatible with `numpy.asanyarray` is possible
    (generally producing ``NdarrayDistribution``).

    See also: https://docs.astropy.org/en/stable/uncertainty/

    Parameters
    ----------
    samples : array-like
        The distribution, with sampling along the *trailing* axis. If 1D, the sole
        dimension is used as the sampling axis (i.e., it is a scalar distribution).
        If an |ndarray| or subclass, the data will not be copied unless it is not
        possible to take a view (generally, only when the strides of the last axis
        are negative).

    """

    _generated_subclasses = {}

    def __new__(cls, samples):
        if isinstance(samples, Distribution):
            samples = samples.distribution
        else:
            # We can handle strides larger than the itemsize by stride trickery,
            # but the samples do have to have positive offset from each other,
            # so if that's not the case, we make a copy.
            if not isinstance(samples, np.ndarray) or (
                samples.strides[-1] < samples.dtype.itemsize
                and samples.shape[-1:] != (1,)  # strides[-1] can be 0.
            ):
                samples = np.asanyarray(samples, order="C")
            if samples.shape == ():
                raise TypeError("Attempted to initialize a Distribution with a scalar")

        # Do the view in two stages, since for viewing as a new class, it is good
        # to have with the dtype in place (e.g., for Quantity.__array_finalize__,
        # it is needed to create the correct StructuredUnit).
        # First, get the structured dtype with a single "samples" field that holds
        # a size n_samples array of "sample" fields with individual samples.
        # We need the double indirection to be able to deal with non-contiguous samples
        # (such as result from indexing a structured array).
        # Note that for a single sample, samples.strides[-1] can be 0.
        new_dtype = cls._get_distribution_dtype(
            samples.dtype, samples.shape[-1], itemsize=samples.strides[-1]
        )
        # Create a structured array with the new dtype. We ensure we start with a
        # regular ndarray, to avoid interference between our complicated dtype and
        # something else (such as would be the case for Quantity.__array_finalize__,
        # which needs the dtype to be correct for creating a StructuredUnit).
        if samples.dtype.itemsize == samples.strides[-1]:
            # For contiguous last axis, a plain new view suffices.
            structured = samples.view(np.ndarray)
        else:
            # But if the last axis is not contiguous, we use __array_interface__ (as in
            # np.lib.stridetricks.as_strided) to create a new ndarray where the last
            # axis is covered by a void of the same size as the structured new_dtype
            # (assigning it directly does not work: possible empty fields get entries).
            interface = dict(samples.__array_interface__)
            interface["descr"] = interface["typestr"] = f"|V{new_dtype.itemsize}"
            interface["shape"] = samples.shape[:-1]
            interface["strides"] = samples.strides[:-1]
            structured = np.asarray(DummyArray(interface, base=samples))
        # Set our new structured dtype.
        structured.dtype = new_dtype
        # Get rid of trailing dimension of 1.
        structured.shape = samples.shape[:-1]

        # Now view as the Distribution subclass, and finalize based on the
        # original samples (e.g., to set the unit for QuantityDistribution).
        new_cls = cls._get_distribution_cls(type(samples))
        self = structured.view(new_cls)
        if not NUMPY_LT_1_23 or callable(self.__array_finalize__):
            self.__array_finalize__(samples)
        return self

    @classmethod
    def _get_distribution_cls(cls, samples_cls):
        if issubclass(samples_cls, Distribution):
            return samples_cls

        new_cls = cls._generated_subclasses.get(samples_cls)
        if new_cls is None:
            # Walk through MRO and find closest base data class.
            # Note: right now, will basically always be ndarray, but
            # one could imagine needing some special care for one subclass,
            # which would then get its own entry.  E.g., if MaskedAngle
            # defined something special, then MaskedLongitude should depend
            # on it.
            for mro_item in samples_cls.__mro__:
                base_cls = cls._generated_subclasses.get(mro_item)
                if base_cls is not None:
                    break
            else:
                # Just hope that NdarrayDistribution can handle it.
                # TODO: this covers the case where a user puts in a list or so,
                # but for those one could just explicitly do something like
                # _generated_subclasses[list] = NdarrayDistribution
                return NdarrayDistribution

            # We need to get rid of the top _DistributionRepr, since we add
            # it again below (it should always be on top).
            # TODO: remove the need for _DistributionRepr by defining
            # __array_function__ and overriding array2string.
            if base_cls.__mro__[1] is not _DistributionRepr:
                # Sanity check. This should never happen!
                raise RuntimeError(
                    f"found {base_cls=}, which does not have _DistributionRepr at "
                    "__mro__[1]. Please raise an issue at "
                    "https://github.com/astropy/astropy/issues/new/choose."
                )
            # Create (and therefore register) new Distribution subclass for the
            # given samples_cls.
            new_cls = type(
                samples_cls.__name__ + "Distribution",
                (_DistributionRepr, samples_cls) + base_cls.__mro__[2:],
                {"_samples_cls": samples_cls},
            )
            cls._generated_subclasses[samples_cls] = new_cls

        return new_cls

    @classmethod
    def _get_distribution_dtype(cls, dtype, n_samples, itemsize=None):
        dtype = np.dtype(dtype)
        # If not a sample dtype already, create one with a "samples" entry that holds
        # all the samples.  Here, we create an indirection via "sample" because
        # selecting an item from a structured array will necessarily give
        # non-consecutive samples, and these can only be dealt with samples that are
        # larger than a single element (plus stride tricks; see __new__). For those
        # cases, itemsize larger than dtype.itemsize will be passed in (note that for
        # the n_sample=1 case, itemsize can be 0; like itemsize=None, this will be dealt
        # with by "or dtype.itemsize" below).
        if dtype.names != ("samples",):
            sample_dtype = np.dtype(
                dict(
                    names=["sample"],
                    formats=[dtype],
                    itemsize=itemsize or dtype.itemsize,
                )
            )
            dtype = np.dtype([("samples", sample_dtype, (n_samples,))])
        return dtype

    @property
    def distribution(self):
        return self["samples"]["sample"]

    @property
    def dtype(self):
        return super().dtype["samples"].base["sample"]

    @dtype.setter
    def dtype(self, dtype):
        dtype = self._get_distribution_dtype(
            dtype, self.n_samples, itemsize=super().dtype["samples"].base.itemsize
        )
        super(Distribution, self.__class__).dtype.__set__(self, dtype)

    def astype(self, dtype, *args, **kwargs):
        dtype = self._get_distribution_dtype(dtype, self.n_samples)
        return super().astype(dtype, *args, **kwargs)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        converted = []
        outputs = kwargs.pop("out", None)
        if outputs:
            kwargs["out"] = tuple(
                (output.distribution if isinstance(output, Distribution) else output)
                for output in outputs
            )
            if ufunc.nout == 1:
                outputs = outputs[0]

        axis = kwargs.get("axis", None)
        keepdims = kwargs.get("keepdims", False)
        if method in {"reduce", "accumulate", "reduceat"}:
            if axis is None:
                assert isinstance(inputs[0], Distribution)
                kwargs["axis"] = tuple(range(inputs[0].ndim))

        for input_ in inputs:
            # For any input that is not a Distribution, we add an axis at the
            # end, to allow proper broadcasting with the distributions.
            if isinstance(input_, Distribution):
                converted.append(input_.distribution)
            else:
                shape = getattr(input_, "shape", ())
                if shape:
                    converted.append(input_[..., np.newaxis])
                else:
                    converted.append(input_)

        if ufunc.signature:
            # We're dealing with a gufunc, which may add an extra axis.
            # Ignore axes keyword for now...  TODO: remove this limitation.
            if "axes" in kwargs:
                raise NotImplementedError(
                    "Distribution does not yet support gufunc calls with 'axes'."
                )
            # Parse signature with private numpy function. Note it
            # cannot handle spaces in tuples, so remove those.
            in_sig, out_sig = _parse_gufunc_signature(ufunc.signature.replace(" ", ""))
            ncore_in = [len(sig) for sig in in_sig]
            ncore_out = [len(sig) + (1 if keepdims else 0) for sig in out_sig]
            if ufunc.nout == 1:
                ncore_out = ncore_out[0]
            if axis is None:
                axis = -1
            converted = [
                np.moveaxis(conv, -1, normalize_axis_index(axis, conv.ndim) - ncore)
                if ncore and getattr(conv, "shape", ())
                else conv
                for conv, ncore in zip(converted, ncore_in)
            ]
        else:
            ncore_out = None

        results = getattr(ufunc, method)(*converted, **kwargs)

        return self._result_as_distribution(results, outputs, ncore_out, axis)

    def __array_function__(self, function, types, args, kwargs):
        # TODO: go through functions systematically to see which ones
        # work and/or can be supported.
        if function in FUNCTION_HELPERS:
            function_helper = FUNCTION_HELPERS[function]
            try:
                args, kwargs, out = function_helper(*args, **kwargs)
            except NotImplementedError:
                return self._not_implemented_or_raise(function, types)

            result = super().__array_function__(function, types, args, kwargs)
            # Fall through to return section, wrapping distributions, unless
            # indicated otherwise.  TODO: use a less hacky solution?
            if out is True:
                return result

        else:  # pragma: no cover
            # By default, just pass it through for now.
            return super().__array_function__(function, types, args, kwargs)

        # We're done if the result was NotImplemented, which can happen
        # if other inputs/outputs override __array_function__;
        # hopefully, they can then deal with us.
        if result is NotImplemented:
            return NotImplemented

        return self._result_as_distribution(result, out=out)

    def _result_as_distribution(self, result, out, ncore_out=None, axis=None):
        """Turn result into a distribution.

        If no output is given, it will create a Distribution from the array,
        If an output is given, it should be fine as is.

        Parameters
        ----------
        result : ndarray or tuple thereof
            Array(s) which need to be turned into Distribution.
        out : Distribution, tuple of Distribution or None
            Possible output |Distribution|. Should be `None` or a tuple if result
            is a tuple.
        ncore_out: int or tuple thereof
            The number of core dimensions for the output array for a gufunc.  This
            is used to determine which axis should be used for the samples.
        axis: int or None
            The axis a gufunc operated on.  Used only if ``ncore_out`` is given.

        Returns
        -------
        out : Distribution
        """
        if isinstance(result, (tuple, list)):
            if out is None:
                out = (None,) * len(result)
            if ncore_out is None:
                ncore_out = (None,) * len(result)
            # Some np.linalg functions return namedtuple, which is handy to access
            # elements by name, but cannot be directly initialized with an iterator.
            result_cls = getattr(result, "_make", result.__class__)
            return result_cls(
                self._result_as_distribution(result_, out_, axis=axis, ncore_out=ncore)
                for (result_, out_, ncore) in zip(result, out, ncore_out)
            )

        if out is None:
            # Turn the result into a Distribution if needed.
            if not isinstance(result, Distribution) and getattr(result, "shape", ()):
                if ncore_out is not None:
                    result = np.moveaxis(
                        result, normalize_axis_index(axis, result.ndim) - ncore_out, -1
                    )
                return Distribution(result)
            else:
                return result
        else:
            # TODO: remove this sanity check once test cases are more complete.
            assert isinstance(out, Distribution)
            return out

    def _not_implemented_or_raise(self, function, types):
        # Our function helper or dispatcher found that the function does not
        # work with Distribution.  In principle, there may be another class that
        # knows what to do with us, for which we should return NotImplemented.
        # But if there is ndarray (or a non-Distribution subclass of it) around,
        # it quite likely coerces, so we should just break.
        if any(
            issubclass(t, np.ndarray) and not issubclass(t, Distribution) for t in types
        ):
            raise TypeError(
                f"the Distribution implementation cannot handle {function} "
                "with the given arguments."
            ) from None
        else:
            return NotImplemented

    # Override __eq__ and __ne__ to pass on directly to the ufunc since
    # otherwise comparisons with non-distributions do not work (but
    # deferring if other defines __array_ufunc__ = None -- see
    # numpy/core/src/common/binop_override.h for the logic; we assume we
    # will never deal with __array_priority__ any more).  Note: there is no
    # problem for other comparisons, since for those, structured arrays are
    # not treated differently in numpy/core/src/multiarray/arrayobject.c.
    def __eq__(self, other):
        if getattr(other, "__array_ufunc__", False) is None:
            return NotImplemented
        return np.equal(self, other)

    def __ne__(self, other):
        if getattr(other, "__array_ufunc__", False) is None:
            return NotImplemented
        return np.not_equal(self, other)

    @property
    def n_samples(self):
        """
        The number of samples of this distribution.  A single `int`.
        """
        return super().dtype["samples"].shape[0]

    def pdf_mean(self, dtype=None, out=None):
        """
        The mean of this distribution.

        Arguments are as for `numpy.mean`.
        """
        return self.distribution.mean(axis=-1, dtype=dtype, out=out)

    def pdf_std(self, dtype=None, out=None, ddof=0):
        """
        The standard deviation of this distribution.

        Arguments are as for `numpy.std`.
        """
        return self.distribution.std(axis=-1, dtype=dtype, out=out, ddof=ddof)

    def pdf_var(self, dtype=None, out=None, ddof=0):
        """
        The variance of this distribution.

        Arguments are as for `numpy.var`.
        """
        return self.distribution.var(axis=-1, dtype=dtype, out=out, ddof=ddof)

    def pdf_median(self, out=None):
        """
        The median of this distribution.

        Parameters
        ----------
        out : array, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output,
            but the type (of the output) will be cast if necessary.
        """
        return np.median(self.distribution, axis=-1, out=out)

    def pdf_mad(self, out=None):
        """
        The median absolute deviation of this distribution.

        Parameters
        ----------
        out : array, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output,
            but the type (of the output) will be cast if necessary.
        """
        median = self.pdf_median(out=out)
        absdiff = np.abs(self - median)
        return np.median(
            absdiff.distribution, axis=-1, out=median, overwrite_input=True
        )

    def pdf_smad(self, out=None):
        """
        The median absolute deviation of this distribution rescaled to match the
        standard deviation for a normal distribution.

        Parameters
        ----------
        out : array, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output,
            but the type (of the output) will be cast if necessary.
        """
        result = self.pdf_mad(out=out)
        result *= SMAD_SCALE_FACTOR
        return result

    def pdf_percentiles(self, percentile, **kwargs):
        """
        Compute percentiles of this Distribution.

        Parameters
        ----------
        percentile : float or array of float or `~astropy.units.Quantity`
            The desired percentiles of the distribution (i.e., on [0,100]).
            `~astropy.units.Quantity` will be converted to percent, meaning
            that a ``dimensionless_unscaled`` `~astropy.units.Quantity` will
            be interpreted as a quantile.

        Additional keywords are passed into `numpy.percentile`.

        Returns
        -------
        percentiles : `~astropy.units.Quantity` ['dimensionless']
            The ``fracs`` percentiles of this distribution.
        """
        percentile = u.Quantity(percentile, u.percent).value
        percs = np.percentile(self.distribution, percentile, axis=-1, **kwargs)
        # numpy.percentile strips units for unclear reasons, so we have to make
        # a new object with units
        if hasattr(self.distribution, "_new_view"):
            return self.distribution._new_view(percs)
        else:
            return percs

    def pdf_histogram(self, **kwargs):
        """
        Compute histogram over the samples in the distribution.

        Parameters
        ----------
        All keyword arguments are passed into `astropy.stats.histogram`. Note
        That some of these options may not be valid for some multidimensional
        distributions.

        Returns
        -------
        hist : array
            The values of the histogram. Trailing dimension is the histogram
            dimension.
        bin_edges : array of dtype float
            Return the bin edges ``(length(hist)+1)``. Trailing dimension is the
            bin histogram dimension.
        """
        distr = self.distribution
        raveled_distr = distr.reshape(distr.size // distr.shape[-1], distr.shape[-1])

        nhists = []
        bin_edges = []
        for d in raveled_distr:
            nhist, bin_edge = stats.histogram(d, **kwargs)
            nhists.append(nhist)
            bin_edges.append(bin_edge)

        nhists = np.array(nhists)
        nh_shape = self.shape + (nhists.size // self.size,)
        bin_edges = np.array(bin_edges)
        be_shape = self.shape + (bin_edges.size // self.size,)
        return nhists.reshape(nh_shape), bin_edges.reshape(be_shape)


class ScalarDistribution(Distribution, np.void):
    """Scalar distribution.

    This class mostly exists to make `~numpy.array2print` possible for
    all subclasses.  It is a scalar element, still with n_samples samples.
    """

    pass


class ArrayDistribution(Distribution, np.ndarray):
    # This includes the important override of view and __getitem__
    # which are needed for all ndarray subclass Distributions, but not
    # for the scalar one.
    _samples_cls = np.ndarray

    # Override view so that we stay a Distribution version of the new type.
    def view(self, dtype=None, type=None):
        """New view of array with the same data.

        Like `~numpy.ndarray.view` except that the result will always be a new
        `~astropy.uncertainty.Distribution` instance.  If the requested
        ``type`` is a `~astropy.uncertainty.Distribution`, then no change in
        ``dtype`` is allowed.

        """
        if type is None:
            if isinstance(dtype, builtins.type) and issubclass(dtype, np.ndarray):
                type = self._get_distribution_cls(dtype)
                dtype = None
            else:
                type = self.__class__
        else:
            type = self._get_distribution_cls(type)

        type = self._get_distribution_cls(type)
        if dtype is None:
            return super().view(type=type)

        dtype = np.dtype(dtype)
        if dtype == self.dtype:
            return super().view(type=type)

        if dtype.names == ("samples",):
            # Assume the user knows what they are doing.
            return super().view(dtype, type)

        if dtype.shape == () and dtype.itemsize == self.dtype.itemsize:
            dtype = self._get_distribution_dtype(
                dtype,
                self.n_samples,
                itemsize=super(Distribution, self).dtype["samples"].base.itemsize,
            )
            return super().view(dtype, type)

        samples_cls = type._samples_cls
        if dtype.itemsize == self.dtype.itemsize:
            distr = self.distribution
            distr_view = distr.view(dtype, samples_cls)
            # This does not necessarily leave the sample axis in the right place.
            return Distribution(np.moveaxis(distr_view, distr.ndim - 1, -1))
        elif dtype.itemsize == super(Distribution, self).dtype["samples"].base.itemsize:
            distr = np.moveaxis(self.distribution, -1, -2)
            distr_view = distr.view(dtype, samples_cls).squeeze(-1)
            return Distribution(distr_view)
        else:
            raise ValueError(
                f"{self.__class__} can only be viewed with a dtype with "
                "itemsize {self.strides[-1]} or {self.dtype.itemsize}"
            )

    @property
    def distribution(self):
        # Like in the creation, we go through an ndarray to ensure we have our
        # actual dtype and to avoid entering, e.g., Quantity.__getitem__, which
        # would give problems with units.
        structured = super().view(np.ndarray)
        distribution = structured["samples"]["sample"].view(self._samples_cls)
        if not NUMPY_LT_1_23 or callable(self.__array_finalize__):
            distribution.__array_finalize__(self)
        return distribution

    def __getitem__(self, item):
        if isinstance(item, str):
            # "samples" should always get back to the samples class.
            if item == "samples":
                return self.distribution
            else:
                # Hard to get this right directly, so instead get item from the
                # distribution, and create a new instance.  We move the sample axis to
                # the end to ensure the order is right for possible subarrays.
                return Distribution(np.moveaxis(self.distribution[item], self.ndim, -1))

        if isinstance(item, Distribution):
            # Required for in-place operations like dist[dist < 0] += 360.
            return self.distribution[item.distribution]

        result = super().__getitem__(item)
        if isinstance(result, np.void):
            return result.view((ScalarDistribution, result.dtype))
        else:
            return result

    def __setitem__(self, item, value):
        if isinstance(item, Distribution):
            # Support operations like dist[dist < 0] = 0.
            self.distribution[item.distribution] = value
            return

        if isinstance(item, str):
            if item == "samples":
                self.distribution[()] = value
                return

            # Get a view of this item (non-trivial; see above).
            self = self[item]
            item = ()

        if not isinstance(value, Distribution):
            # If value is not already a Distribution, first make it an array
            # to help interpret possible structured dtype, and then turn it
            # into a Distribution with n_samples=1 (which will broadcast).
            value = np.asanyarray(value, dtype=self.dtype)
            value = Distribution(value[..., np.newaxis])

        super().__setitem__(item, value)


class _DistributionRepr:
    def __repr__(self):
        reprarr = repr(self.distribution)
        if reprarr.endswith(">"):
            firstspace = reprarr.find(" ")
            reprarr = reprarr[firstspace + 1 : -1]  # :-1] removes the ending '>'
            return (
                f"<{self.__class__.__name__} {reprarr} with n_samples={self.n_samples}>"
            )
        else:  # numpy array-like
            firstparen = reprarr.find("(")
            reprarr = reprarr[firstparen:]
            return f"{self.__class__.__name__}{reprarr} with n_samples={self.n_samples}"

    def __str__(self):
        distrstr = str(self.distribution)
        toadd = f" with n_samples={self.n_samples}"
        return distrstr + toadd

    def _repr_latex_(self):
        if hasattr(self.distribution, "_repr_latex_"):
            superlatex = self.distribution._repr_latex_()
            toadd = rf", \; n_{{\rm samp}}={self.n_samples}"
            return superlatex[:-1] + toadd + superlatex[-1]
        else:
            return None


class NdarrayDistribution(_DistributionRepr, ArrayDistribution):
    pass


# Ensure our base NdarrayDistribution is known.
Distribution._generated_subclasses[np.ndarray] = NdarrayDistribution
