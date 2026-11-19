# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import warnings

import numpy as np

from astropy.io import fits
from astropy.utils.exceptions import (
    AstropyDeprecationWarning,
    AstropyUserWarning,
)

from ._wcs import DistortionLookupTable

# Regular expression matching distortion paper (WCS Distortion Paper Draft)
# prior-distortion keywords: ``CPDISja``, ``CPERRja`` and record-valued
# ``DPja.<field>``. Groups are the keyword family, the axis number, the
# alternate WCS key ('' for the primary WCS) and the record field.
_PRIOR_DIS_KW = re.compile(r"(CPDIS|CPERR|DP)([0-9]+)([A-Z]?)(?:\.(.+))?")

# HST/STScI convention for prior-distortion lookup tables shared between the
# primary and alternate WCSs (see the ``shared_cpdis`` parameter).
_STSCI_TELESCOPE = "HST"
_STSCI_INSTRUMENTS = frozenset({"ACS", "WFC3"})
_STSCI_AXES = frozenset({1, 2})
_STSCI_MAX_EXTVER = 4  # 2 chips x 2 axes (ACS/WFC, WFC3/UVIS)


class _LookupTableMixin:
    """
    Reading and writing of the lookup-table distortions of a WCS
    (``CPDISja`` and ``D2IMDISj`` keywords, ``WCSDVARR`` and ``D2IMARR``
    extensions).

    Mixed into `~astropy.wcs.WCS`. Uses ``naxis``, ``wcs``, ``cpdis1``,
    ``cpdis2``, ``det2im1`` and ``det2im2``, which ``WCSBase`` provides.
    """

    def _read_det2im_kw(self, header, fobj, err=0.0):
        """
        Create a `distortion paper`_ type lookup table for detector to
        image plane correction.
        """
        if fobj is None:
            return (None, None)

        if not isinstance(fobj, fits.HDUList):
            return (None, None)

        try:
            axiscorr = header["AXISCORR"]
            d2imdis = self._read_d2im_old_format(header, fobj, axiscorr)
            return d2imdis
        except KeyError:
            pass

        dist = "D2IMDIS"
        d_kw = "D2IM"
        err_kw = "D2IMERR"
        tables = {}
        for i in range(1, self.naxis + 1):
            d_error = header.get(err_kw + str(i), 0.0)
            if d_error < err:
                # Skip this table, and remove its keywords as when it is read.
                tables[i] = None
                stale = [
                    kw
                    for kw in header
                    if kw == dist + str(i) or kw.startswith(f"{d_kw}{i}.")
                ]
                for kw in stale:
                    header.remove(kw, ignore_missing=True, remove_all=True)
                continue
            distortion = dist + str(i)
            if distortion in header:
                dis = header[distortion].lower()
                if dis == "lookup":
                    del header[distortion]
                    assert isinstance(fobj, fits.HDUList), (
                        "An astropy.io.fits.HDUList"
                        "is required for Lookup table distortion."
                    )
                    dp = (d_kw + str(i)).strip()
                    dp_extver_key = dp + ".EXTVER"
                    if dp_extver_key in header:
                        d_extver = header[dp_extver_key]
                        del header[dp_extver_key]
                    else:
                        d_extver = 1
                    dp_axis_key = dp + f".AXIS.{i:d}"
                    if i == header[dp_axis_key]:
                        d_data = fobj["D2IMARR", d_extver].data
                    else:
                        d_data = (fobj["D2IMARR", d_extver].data).transpose()
                    del header[dp_axis_key]
                    d_header = fobj["D2IMARR", d_extver].header
                    d_crpix = (d_header.get("CRPIX1", 0.0), d_header.get("CRPIX2", 0.0))
                    d_crval = (d_header.get("CRVAL1", 0.0), d_header.get("CRVAL2", 0.0))
                    d_cdelt = (d_header.get("CDELT1", 1.0), d_header.get("CDELT2", 1.0))
                    d_lookup = DistortionLookupTable(d_data, d_crpix, d_crval, d_cdelt)
                    tables[i] = d_lookup
                else:
                    warnings.warn(
                        "Polynomial distortion is not implemented.\n",
                        AstropyUserWarning,
                    )
                for key in set(header):
                    if key.startswith(dp + "."):
                        del header[key]
            else:
                tables[i] = None
        if not tables:
            return (None, None)
        else:
            return (tables.get(1), tables.get(2))

    def _read_d2im_old_format(self, header, fobj, axiscorr):
        warnings.warn(
            "The use of ``AXISCORR`` for D2IM correction has been"
            " deprecated.`~astropy.wcs` will read in files with ``AXISCORR`` but"
            " ``to_fits()`` will write out files without it.",
            AstropyDeprecationWarning,
        )
        cpdis = [None, None]
        crpix = [0.0, 0.0]
        crval = [0.0, 0.0]
        cdelt = [1.0, 1.0]
        try:
            d2im_data = fobj[("D2IMARR", 1)].data
        except KeyError:
            return (None, None)
        except AttributeError:
            return (None, None)

        d2im_data = np.array([d2im_data])
        d2im_hdr = fobj[("D2IMARR", 1)].header
        naxis = d2im_hdr["NAXIS"]

        for i in range(1, naxis + 1):
            crpix[i - 1] = d2im_hdr.get("CRPIX" + str(i), 0.0)
            crval[i - 1] = d2im_hdr.get("CRVAL" + str(i), 0.0)
            cdelt[i - 1] = d2im_hdr.get("CDELT" + str(i), 1.0)

        cpdis = DistortionLookupTable(d2im_data, crpix, crval, cdelt)

        if axiscorr == 1:
            return (cpdis, None)
        elif axiscorr == 2:
            return (None, cpdis)
        else:
            warnings.warn("Expected AXISCORR to be 1 or 2", AstropyUserWarning)
            return (None, None)

    def _write_det2im(self, hdulist, ext=0):
        """
        Write the detector-to-image lookup tables (``D2IMDISj``) to
        ``hdulist[ext]``. ``D2IMDISj`` keywords have no alternate WCS key:
        they are shared by all WCSs in a header.
        """
        _write_lookup_tables(
            hdulist,
            ext,
            (self.det2im1, self.det2im2),
            dist="D2IMDIS",
            rec_kw="D2IM",
            err_kw="D2IMERR",
            extname="D2IMARR",
            label="D2IM",
            type_comment="Detector to image correction type",
            key="",
        )
        if self.det2im1 is not None or self.det2im2 is not None:
            # A header in the old format (AXISCORR) is read in that format
            # even if it has D2IMDISj keywords, so remove the old keyword.
            hdulist[ext].header.remove("AXISCORR", ignore_missing=True)

    def _read_distortion_kw(
        self, header, fobj, *, dist="CPDIS", err=0.0, key=" ", inherit_primary=False
    ):
        """
        Reads `distortion paper`_ table-lookup keywords and data, and
        returns a 2-tuple of `~astropy.wcs.DistortionLookupTable`
        objects.
        If no `distortion paper`_ keywords are found, ``(None, None)``
        is returned.

        Parameters
        ----------
        header : `~astropy.io.fits.Header` or str or bytes
            FITS header containing distortion keywords
        fobj : `~astropy.io.fits.HDUList` or None
            FITS file object for lookup table data
        dist : str, optional
            Distortion type: "CPDIS" for prior or "CQDIS" for sequent.
            Default is "CPDIS".
        err : float, optional
            Minimum error threshold. Distortions with maximum error less
            than this value will be ignored. Default is 0.0.
        key : str, optional
            WCS key character (e.g., ' ', 'A', 'B', etc.). This is used
            to construct the proper keyword names with the appropriate
            suffix. Default is ' ' (primary WCS).
        inherit_primary : bool, optional
            If `True`, read the unsuffixed (primary WCS) keywords instead of
            those suffixed with *key*. Default is `False`.

        Returns
        -------
        tuple
            A 2-tuple of `~astropy.wcs.DistortionLookupTable` objects
            for axes 1 and 2, or ``(None, None)`` if no distortion
            keywords are found.
        """
        if isinstance(header, (str, bytes)):
            return (None, None)

        key_suffix = "" if inherit_primary else key.strip().upper()

        if dist == "CPDIS":
            d_kw = "DP"
            err_kw = "CPERR"
        else:
            d_kw = "DQ"
            err_kw = "CQERR"

        tables = {}
        for i in range(1, self.naxis + 1):
            axis_key_suffix = str(i) + key_suffix
            d_error_key = err_kw + axis_key_suffix

            if d_error_key in header:
                d_error = header[d_error_key]
                del header[d_error_key]
            else:
                d_error = 0.0

            if d_error < err:
                # Skip this table, but remove its keywords: WCSLIB cannot
                # evaluate a Lookup distortion if they reach it.
                tables[i] = None
                dis_kw = dist + axis_key_suffix
                dp_prefix = d_kw + axis_key_suffix + "."
                stale = [
                    kw for kw in header if kw == dis_kw or kw.startswith(dp_prefix)
                ]
                for kw in stale:
                    header.remove(kw, ignore_missing=True, remove_all=True)
                continue

            distortion = dist + axis_key_suffix

            if distortion in header:
                dis = header[distortion].lower()
                del header[distortion]

                if dis == "lookup":
                    if not isinstance(fobj, fits.HDUList):
                        raise ValueError(
                            "an astropy.io.fits.HDUList is "
                            "required for Lookup table distortion."
                        )

                    dp = d_kw + axis_key_suffix

                    dp_extver_key = dp + ".EXTVER"
                    if dp_extver_key in header:
                        d_extver = header[dp_extver_key]
                        del header[dp_extver_key]
                    else:
                        d_extver = 1

                    dp_axis_key = dp + f".AXIS.{i:d}"
                    if dp_axis_key in header:
                        axis_val = header[dp_axis_key]
                        if i == axis_val:
                            d_data = fobj["WCSDVARR", d_extver].data
                        else:
                            d_data = (fobj["WCSDVARR", d_extver].data).transpose()
                        del header[dp_axis_key]
                    else:
                        d_data = fobj["WCSDVARR", d_extver].data

                    d_header = fobj["WCSDVARR", d_extver].header
                    d_crpix = (d_header.get("CRPIX1", 0.0), d_header.get("CRPIX2", 0.0))
                    d_crval = (d_header.get("CRVAL1", 0.0), d_header.get("CRVAL2", 0.0))
                    d_cdelt = (d_header.get("CDELT1", 1.0), d_header.get("CDELT2", 1.0))
                    d_lookup = DistortionLookupTable(d_data, d_crpix, d_crval, d_cdelt)
                    tables[i] = d_lookup

                    for kwd in set(header):
                        if kwd.startswith(dp + "."):
                            del header[kwd]
                else:
                    warnings.warn(
                        "Polynomial distortion is not implemented.\n",
                        AstropyUserWarning,
                    )
            else:
                tables[i] = None

        if not tables:
            return (None, None)
        else:
            return (tables.get(1), tables.get(2))

    def _write_distortion_kw(self, hdulist, dist="CPDIS", ext=0, key=" "):
        """
        Write the `distortion paper`_ lookup tables (``CPDISja``, or
        ``CQDISia`` for ``dist="CQDIS"``) for WCS ``key`` to ``hdulist[ext]``.
        """
        if dist == "CPDIS":
            rec_kw = "DP"
            err_kw = "CPERR"
            comment = "Prior distortion function type"
        else:
            rec_kw = "DQ"
            err_kw = "CQERR"
            comment = "Sequent distortion function type"
        _write_lookup_tables(
            hdulist,
            ext,
            (self.cpdis1, self.cpdis2),
            dist=dist,
            rec_kw=rec_kw,
            err_kw=err_kw,
            extname="WCSDVARR",
            label=dist,
            type_comment=comment,
            key=key,
        )

    def write_lookup_tables(self, hdulist, ext=0, key=None):
        """
        Write the lookup-table distortions of this WCS to an existing
        `~astropy.io.fits.HDUList`.

        The `distortion paper`_ prior distortion (`cpdis1`, `cpdis2`) is
        written as ``CPDISja``, ``CPERRja`` and ``DPja`` keywords for WCS
        *key*, and the detector-to-image correction (`det2im1`, `det2im2`) as
        ``D2IMDISj``, ``D2IMERRj`` and ``D2IMj`` keywords, to the header of
        ``hdulist[ext]``. The tables themselves are stored in ``WCSDVARR``
        and ``D2IMARR`` image extensions of *hdulist*.

        Use it together with `to_header` to add this WCS to a file that
        already holds other WCSs, e.g.::

            hdr = hdulist[ext].header
            hdr.update(w.to_header(key="A", relax=True))
            w.write_lookup_tables(hdulist, ext, key="A")

        Parameters
        ----------
        hdulist : `~astropy.io.fits.HDUList`
            The file to write to. Image extensions are appended to it as
            needed.
        ext : int, str or tuple, optional
            The HDU whose header receives the keywords. Default is 0.
        key : str, optional
            The WCS key, ``' '`` or ``'A'``-``'Z'``. Defaults to the key of
            this WCS (``self.wcs.alt``).

        Notes
        -----
        - The lookup-table keywords (``CPDISja``, ``CPERRja`` and ``DPja``)
          that the header already has for *key* are replaced, and removed
          for an axis that has no table in this WCS. Keywords of other keys
          are not changed. Neither are keywords of *key* that describe
          another type of prior distortion, on an axis that gets no table.
          A header that already describes the tables of this WCS for *key*
          is left as it is.

        - Each table is stored in a ``WCSDVARR`` (or ``D2IMARR``) extension
          of its own. A table that replaces another table of the same *key*
          and axis is written over it, with the same ``EXTVER``, unless
          another keyword in *hdulist* refers to that extension. Otherwise,
          a table that is identical to an extension already in *hdulist* is
          not written again; the keywords refer to the existing extension.
          Any other table is appended, with ``EXTVER`` equal to the axis
          number if it is free and the smallest free ``EXTVER`` otherwise.
          Extensions are never removed.

        - ``D2IMDISj`` keywords have no WCS key and apply to all WCSs in the
          header. This WCS's detector-to-image correction is written if the
          header has none. If the header already has one that differs, it is
          replaced when *key* is ``' '`` and a `ValueError` is raised
          otherwise. If this WCS has no detector-to-image correction,
          existing ``D2IMDISj`` keywords are left unchanged. An ``AXISCORR``
          keyword (the old format of this correction) is removed when the
          correction is written.

        - ``CPERRja`` and ``D2IMERRj`` are set to the largest absolute value
          in each table (ignoring NaN and infinite values), so that the
          ``minerr`` argument of `WCS` works when the file is read back.

        - In HST ACS and WFC3 files that follow the STScI convention of
          sharing the unsuffixed prior-distortion tables with all alternate
          WCSs (see the ``shared_cpdis`` parameter of `WCS`), writing keyed
          ``CPDISja`` keywords turns that sharing off for the whole header,
          and changing the unsuffixed tables (*key* ``' '``) changes them for
          every alternate WCS that shares them. An
          `~astropy.utils.exceptions.AstropyUserWarning` is issued in both
          cases. It also says so if the header no longer follows the
          convention after the unsuffixed tables are changed, for example
          because the new tables could not be written over the old ones and
          no ``EXTVER`` from 1 to 4 was free.
        """
        if key is None:
            key = self.wcs.alt
        key = (key or " ").upper()
        if len(key) != 1 or not (key == " " or "A" <= key <= "Z"):
            raise ValueError("key must be ' ' or 'A'-'Z'")

        header = hdulist[ext].header

        if self.det2im1 is not None or self.det2im2 is not None:
            stored = self._stored_det2im(hdulist, ext)
            has_d2im = stored[0] is not None or stored[1] is not None
            same = _same_lookup_table(stored[0], self.det2im1)
            same = same and _same_lookup_table(stored[1], self.det2im2)
            if not has_d2im or (not same and key == " "):
                self._write_det2im(hdulist, ext)
            elif not same:
                raise ValueError(
                    "The header already has a different detector-to-image "
                    "correction (D2IMDISj). It is shared by all WCSs in the "
                    f"header, so it cannot be changed for WCS key {key!r}; "
                    "write it with key=' ' instead."
                )

        stored = self._stored_cpdis(hdulist, ext, key)
        if (
            stored is not None
            and _same_lookup_table(stored[0], self.cpdis1)
            and _same_lookup_table(stored[1], self.cpdis2)
        ):
            # The header already describes these tables for this key.
            return

        shared = _is_hst_acs_wfc3(hdulist) and _is_stsci_prior_distortion(header)
        self._write_distortion_kw(hdulist, dist="CPDIS", ext=ext, key=key)
        if shared:
            _warn_shared_cpdis_changed(header, key)

    def _stored_det2im(self, hdulist, ext):
        """
        Return the detector-to-image lookup tables that the header of
        ``hdulist[ext]`` describes: ``(None, None)`` if there are none, or
        if they cannot be read.
        """
        header = hdulist[ext].header.copy()
        try:
            with warnings.catch_warnings():
                # Reading the old (AXISCORR) format warns that it is
                # deprecated, which is of no interest to someone writing.
                warnings.simplefilter("ignore", AstropyDeprecationWarning)
                warnings.simplefilter("ignore", AstropyUserWarning)
                return self._read_det2im_kw(header, hdulist)
        except (KeyError, TypeError, ValueError):
            # e.g. the keywords refer to an extension that does not exist
            return (None, None)

    def _stored_cpdis(self, hdulist, ext, key):
        """
        Return the prior-distortion lookup tables that the header of
        ``hdulist[ext]`` has for WCS ``key`` itself (not the tables that an
        alternate WCS shares with the primary WCS), or `None` if they cannot
        be read.
        """
        header = hdulist[ext].header.copy()
        try:
            with warnings.catch_warnings():
                # Another type of distortion only gives a warning when read.
                warnings.simplefilter("ignore", AstropyUserWarning)
                return self._read_distortion_kw(header, hdulist, dist="CPDIS", key=key)
        except (KeyError, TypeError, ValueError):
            # e.g. the keywords refer to an extension that does not exist
            return None


def _inherit_primary_prior_distortion(header, fobj, key, shared_cpdis):
    """
    Decide whether the alternate WCS ``key`` should use the unsuffixed
    (primary WCS) prior-distortion keywords. See the
    ``shared_cpdis`` parameter of `WCS` for details.
    """
    if not isinstance(header, fits.Header):
        return False
    key = (key or " ").strip().upper()
    if not key or (shared_cpdis is not None and not shared_cpdis):
        return False
    if _has_prior_distortion_kw(header, key):
        return False  # the alternate WCS has its own keywords
    if not _has_prior_distortion_kw(header, ""):
        return False  # nothing to share
    if shared_cpdis is None:
        return _is_hst_acs_wfc3(fobj) and _is_stsci_prior_distortion(header)
    return True


def _same_lookup_table(a, b):
    """
    Return `True` if two `~astropy.wcs.DistortionLookupTable` objects (or
    `None`) describe the same correction. NaN values at the same positions
    count as equal.
    """
    if a is None or b is None:
        return a is None and b is None
    return (
        np.array_equal(a.data, b.data, equal_nan=True)
        and np.array_equal(a.crpix, b.crpix)
        and np.array_equal(a.crval, b.crval)
        and np.array_equal(a.cdelt, b.cdelt)
    )


def _as_int(value):
    """Return ``value`` as an `int` if it is a whole number, otherwise `None`."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return int(f) if f.is_integer() else None


def _prior_distortion_keywords(header):
    """
    Yield ``(keyword, family, axis, key, field)`` for every prior-distortion
    keyword in ``header``: the keyword as found in ``header``; ``'CPDIS'``,
    ``'CPERR'`` or ``'DP'``; the axis number; the alternate WCS key (``''``
    for the primary WCS); and the record field of a ``DPja`` keyword (`None`
    otherwise, including for a ``DPja`` card whose value is not a record).
    """
    for keyword in header:
        match = _PRIOR_DIS_KW.fullmatch(str(keyword).strip().upper())
        if match is None:
            continue
        family, axis, key, field = match.groups()
        if field is not None and family != "DP":
            continue  # e.g. "CPDIS1.X" is not a distortion keyword
        yield keyword, family, int(axis), key, field


def _has_prior_distortion_kw(header, key):
    """
    Return `True` if ``header`` contains any ``CPDISja``, ``CPERRja`` or
    ``DPja`` keyword for WCS ``key`` (``''`` for the primary WCS).
    """
    return any(k == key for _, _, _, k, _ in _prior_distortion_keywords(header))


def _is_hst_acs_wfc3(fobj):
    """
    Return `True` if the primary header of ``fobj`` has ``TELESCOP = 'HST'``
    and ``INSTRUME`` equal to ``'ACS'`` or ``'WFC3'``.
    """
    if not isinstance(fobj, fits.HDUList) or not len(fobj):
        return False
    phdr = fobj[0].header
    telescope = str(phdr.get("TELESCOP", "")).strip().upper()
    instrument = str(phdr.get("INSTRUME", "")).strip().upper()
    return telescope == _STSCI_TELESCOPE and instrument in _STSCI_INSTRUMENTS


def _is_stsci_prior_distortion(header):  # noqa: PLR0911
    """
    Return `True` if every prior-distortion keyword in ``header`` follows the
    layout written by STScI for HST non-polynomial (``NPOLFILE``) distortion
    corrections:

    - there are no keyed (alternate WCS) ``CPDISja``, ``CPERRja`` or ``DPja``
      keywords;
    - prior-distortion keywords exist for axes 1 and 2 only, and for both;
    - ``CPDISj`` is ``'Lookup'`` (case-insensitive);
    - ``DPj.AXIS.1`` and ``DPj.AXIS.2`` are present, have the values 1 and 2
      in either order, and there are no other ``DPj.AXIS.n`` fields;
    - ``DPj.EXTVER`` is a whole number from 1 to 4, different for the two
      axes.

    ``CPERRj`` and all other ``DPj`` record fields (``NAXES``, ``OFFSET``,
    ``SCALE``, ...) are ignored.
    """
    axes = {}
    for keyword, family, axis, key, field in _prior_distortion_keywords(header):
        if key or axis not in _STSCI_AXES:
            return False
        entry = axes.setdefault(axis, {"dp": {}})
        if family == "CPDIS":
            entry["cpdis"] = header[keyword]
        elif family == "DP":
            if field is None:
                return False  # a DPj card whose value is not 'field: number'
            entry["dp"][field.upper()] = header[keyword]

    if set(axes) != _STSCI_AXES:
        return False

    extvers = set()
    for entry in axes.values():
        if str(entry.get("cpdis", "")).strip().upper() != "LOOKUP":
            return False
        dp = entry["dp"]
        if {f for f in dp if f.startswith("AXIS.")} != {"AXIS.1", "AXIS.2"}:
            return False
        if {_as_int(dp["AXIS.1"]), _as_int(dp["AXIS.2"])} != _STSCI_AXES:
            return False
        extver = _as_int(dp.get("EXTVER"))
        if extver is None or not 1 <= extver <= _STSCI_MAX_EXTVER:
            return False
        extvers.add(extver)

    return len(extvers) == len(_STSCI_AXES)


def _lookup_table_from_hdu(hdu):
    """
    Lookup table stored in an image extension, as the WCS reader reads it
    for keywords that list the axes of the table in their natural order
    (``AXIS.1`` is 1 and ``AXIS.2`` is 2). `_write_lookup_tables` always
    writes the keywords that way.
    """
    header = hdu.header
    return DistortionLookupTable(
        hdu.data,
        (header.get("CRPIX1", 0.0), header.get("CRPIX2", 0.0)),
        (header.get("CRVAL1", 0.0), header.get("CRVAL2", 0.0)),
        (header.get("CDELT1", 1.0), header.get("CDELT2", 1.0)),
    )


def _hdu_holds_lookup_table(hdu, table):
    """Return `True` if the image extension ``hdu`` holds ``table``."""
    if not isinstance(hdu, fits.ImageHDU):
        return False
    # Compare the shapes first: this does not load the data of the extension.
    if hdu.shape != table.data.shape:
        return False
    if hdu.data is None:
        return False
    return _same_lookup_table(_lookup_table_from_hdu(hdu), table)


def _largest_correction(table):
    """
    Largest absolute value in a lookup table. NaN and infinite values are
    ignored (a FITS header card cannot hold them); 0.0 is returned if the
    table has no finite value.
    """
    magnitudes = np.abs(table.data)
    finite = magnitudes[np.isfinite(magnitudes)]
    if finite.size == 0:
        return 0.0
    return float(finite.max())


def _set_lookup_table_wcs(header, table):
    """Write the reference pixel, value and increment of ``table``."""
    header["CRPIX1"] = (table.crpix[0], "Coordinate system reference pixel")
    header["CRPIX2"] = (table.crpix[1], "Coordinate system reference pixel")
    header["CRVAL1"] = (table.crval[0], "Coordinate system value at reference pixel")
    header["CRVAL2"] = (table.crval[1], "Coordinate system value at reference pixel")
    header["CDELT1"] = (table.cdelt[0], "Coordinate increment along axis")
    header["CDELT2"] = (table.cdelt[1], "Coordinate increment along axis")


def _remove_lookup_keywords(
    header, *, dist, rec_kw, err_kw, suffix, axes_with_new_table=()
):
    """
    Remove from ``header`` the lookup-table keywords (``{dist}ja``,
    ``{err_kw}ja`` and ``{rec_kw}ja``) of WCS key ``suffix``. Return an
    ``{axis: EXTVER}`` dictionary of the lookup tables that the removed
    keywords referred to.

    The keywords of an axis are removed if they describe a lookup table, or
    if the axis is in ``axes_with_new_table`` (a table is about to be
    written for it, whatever the header has now). Keywords of another type
    of distortion (e.g. a ``TPD`` polynomial) are kept on all other axes.
    """
    # Lookup tables that the header has for this key, by axis number.
    old_extvers = {}
    dist_pattern = re.compile(rf"{dist}([0-9]+){suffix}")
    for keyword in header:
        match = dist_pattern.fullmatch(str(keyword).strip().upper())
        if match is None:
            continue
        if str(header[keyword]).strip().upper() != "LOOKUP":
            continue
        axis = int(match.group(1))
        extver = header.get(f"{rec_kw}{axis}{suffix}.EXTVER", 1)
        old_extvers[axis] = _as_int(extver)

    axes_to_clear = set(old_extvers)
    axes_to_clear.update(axes_with_new_table)

    any_pattern = re.compile(rf"(?:{dist}|{err_kw}|{rec_kw})([0-9]+){suffix}(?:\..+)?")
    stale = []
    for index, keyword in enumerate(header):
        match = any_pattern.fullmatch(str(keyword).strip().upper())
        if match is not None and int(match.group(1)) in axes_to_clear:
            stale.append(index)
    for index in reversed(stale):
        del header[index]

    return old_extvers


def _referenced_lookup_extvers(hdulist, *, dist, rec_kw):
    """
    Return the ``EXTVER`` of every lookup-table extension that a header in
    ``hdulist`` refers to: ``{dist}ja = 'LOOKUP'`` with ``{rec_kw}ja.EXTVER``
    (1 if absent), for any axis and WCS key.
    """
    pattern = re.compile(rf"{dist}([0-9]+[A-Z]?)")
    extvers = set()
    for hdu in hdulist:
        header = hdu.header
        if dist == "D2IMDIS" and "AXISCORR" in header:
            # The old format of the detector-to-image correction always
            # uses the first D2IMARR extension.
            extvers.add(1)
        for keyword in header:
            match = pattern.fullmatch(str(keyword).strip().upper())
            if match is None:
                continue
            if str(header[keyword]).strip().upper() != "LOOKUP":
                continue
            extver = header.get(f"{rec_kw}{match.group(1)}.EXTVER", 1)
            extver = _as_int(extver)
            if extver is not None:
                extvers.add(extver)
    return extvers


def _free_extver(preferred, used):
    """
    Return ``preferred`` if it is not in ``used``, and the smallest positive
    integer that is not in ``used`` otherwise.
    """
    if preferred not in used:
        return preferred
    extver = 1
    while extver in used:
        extver += 1
    return extver


def _warn_shared_cpdis_changed(header, key):
    """
    Warn that the prior-distortion keywords just written for ``key``
    change what the alternate WCSs of ``header`` use. ``header`` followed
    the HST/STScI convention of sharing the prior-distortion tables of
    the primary WCS before they were written.
    """
    if key != " ":
        warnings.warn(
            f"Writing CPDISj{key} keywords to a header whose alternate "
            "WCSs share the unsuffixed (CPDISj) prior-distortion tables "
            "of the primary WCS (HST/STScI convention). The other "
            "alternate WCSs will no longer use these tables unless the "
            "file is read with shared_cpdis=True.",
            AstropyUserWarning,
        )
        return

    message = (
        "Changing the unsuffixed (CPDISj) prior-distortion tables of a "
        "header whose alternate WCSs share them (HST/STScI convention). "
    )
    still_shared = _is_stsci_prior_distortion(header)
    tables_removed = not _has_prior_distortion_kw(header, "")
    if still_shared or tables_removed:
        message += (
            "The change applies to all alternate WCSs without CPDISja "
            "keywords of their own."
        )
    else:
        message += (
            "The header no longer follows that convention (tables for "
            "both axes, in WCSDVARR extensions with different EXTVER "
            "from 1 to 4), so the alternate WCSs will not use the new "
            "tables unless the file is read with shared_cpdis=True."
        )
    warnings.warn(message, AstropyUserWarning)


def _write_lookup_tables(
    hdulist,
    ext,
    tables,
    *,
    dist,
    rec_kw,
    err_kw,
    extname,
    label,
    type_comment,
    key,
):
    """
    Write lookup-table distortion keywords for WCS ``key`` to
    ``hdulist[ext].header`` and the tables to ``extname`` image extensions.

    - The lookup-table keywords that the header has for ``key`` are
      replaced, and removed for an axis that has no table in ``tables``;
      see `_remove_lookup_keywords`.
    - ``{err_kw}ja`` is set to the largest absolute value in the table, so
      that the ``minerr`` argument of `WCS` works after a round trip.
    - Each table is stored in an extension of its own, the first of:

      1. the extension that the header already refers to for this key and
         axis, if it holds an identical table;
      2. that same extension, written over with the new table, if no other
         keyword in ``hdulist`` refers to it. Steps 1 and 2 keep the
         ``EXTVER`` that the header has, as the HST/STScI convention for
         sharing tables requires (see `_is_stsci_prior_distortion`);
      3. another extension of ``hdulist`` that holds an identical table;
      4. a new extension, with ``EXTVER`` equal to the axis number if that
         is free (as earlier versions wrote) and the smallest free
         ``EXTVER`` otherwise.

      Two tables written by one call never share an extension, even if they
      are identical: ``WCS.to_fits()`` has always written one extension per
      axis, with ``EXTVER`` equal to the axis number, and code that
      post-processes its output (e.g. ``stwcs``) relies on that.
      Extensions are never removed.
    """
    header = hdulist[ext].header
    suffix = key.strip().upper()

    # Extensions that are in the file before anything is written, by EXTVER.
    # If an EXTVER is repeated, the reader uses the first extension.
    extensions = {}
    for hdu in hdulist:
        if hdu.name == extname and hdu.ver not in extensions:
            extensions[hdu.ver] = hdu

    axes_with_new_table = []
    for num, table in enumerate(tables, start=1):
        if table is not None:
            axes_with_new_table.append(num)

    old_extvers = _remove_lookup_keywords(
        header,
        dist=dist,
        rec_kw=rec_kw,
        err_kw=err_kw,
        suffix=suffix,
        axes_with_new_table=axes_with_new_table,
    )
    # This is done after the keywords of ``key`` are removed, so these are
    # the extensions that something else still uses.
    referenced = _referenced_lookup_extvers(hdulist, dist=dist, rec_kw=rec_kw)

    # EXTVERs given to the tables written by this call.
    taken = set()

    for num, table in enumerate(tables, start=1):
        if table is None:
            continue

        # The extension that the header referred to for this key and axis.
        old_extver = old_extvers.get(num)
        old_hdu = extensions.get(old_extver)
        old_is_available = (
            isinstance(old_hdu, fits.ImageHDU) and old_extver not in taken
        )

        extver = None
        if old_is_available and _hdu_holds_lookup_table(old_hdu, table):
            # 1. The header already refers to this table.
            extver = old_extver
        elif old_is_available and old_extver not in referenced:
            # 2. Nothing else uses the table that is being replaced.
            extver = old_extver
            old_hdu.data = table.data
            _set_lookup_table_wcs(old_hdu.header, table)
        else:
            # 3. Another extension that holds an identical table.
            for candidate, hdu in extensions.items():
                if candidate not in taken and _hdu_holds_lookup_table(hdu, table):
                    extver = candidate
                    break

        if extver is None:
            # 4. A new extension.
            extver = _free_extver(num, set(extensions) | taken)
            image = fits.ImageHDU(table.data, name=extname)
            image.ver = extver
            _set_lookup_table_wcs(image.header, table)
            hdulist.append(image)

        taken.add(extver)

        axis = f"{num:d}{suffix}"
        header[f"{dist}{axis}"] = ("LOOKUP", type_comment)
        header[f"{err_kw}{axis}"] = (
            _largest_correction(table),
            f"Maximum error of {label} correction for axis {num:d}",
        )
        header[f"{rec_kw}{axis}.EXTVER"] = (
            extver,
            f"Version number of {extname} extension",
        )
        header[f"{rec_kw}{axis}.NAXES"] = (
            table.data.ndim,
            f"Number of independent variables in {label} function",
        )
        for i in range(table.data.ndim):
            jth = {1: "1st", 2: "2nd", 3: "3rd"}.get(i + 1, f"{i + 1}th")
            header[f"{rec_kw}{axis}.AXIS.{i + 1:d}"] = (
                i + 1,
                f"Axis number of the {jth} variable in a {label} function",
            )
