/*============================================================================
  WCSLIB 8.3 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2024, Mark Calabretta

  This file is part of WCSLIB.

  WCSLIB is free software: you can redistribute it and/or modify it under the
  terms of the GNU Lesser General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  WCSLIB is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with WCSLIB.  If not, see http://www.gnu.org/licenses.

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: spx.c,v 8.3 2024/05/13 16:33:00 mcalabre Exp $
*===========================================================================*/

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "wcserr.h"
#include "wcsmath.h"
#include "spx.h"


// Map status return value to message.
const char *spx_errmsg[] = {
  "Success",
  "Null spxprm pointer passed",
  "Invalid spectral parameters",
  "Invalid spectral variable",
  "One or more of the inspec coordinates were invalid"};

// Convenience macro for invoking wcserr_set().
#define SPX_ERRMSG(status) WCSERR_SET(status), spx_errmsg[status]

#define C 2.99792458e8
#define h 6.6260755e-34

/*============================================================================
*   Spectral cross conversions; given one spectral coordinate it computes all
*   the others, plus the required derivatives of each with respect to the
*   others.
*===========================================================================*/

int specx(
  const char *type,
  double spec,
  double restfrq,
  double restwav,
  struct spxprm *spx)

{
  static const char *function = "specx";

  if (spx == 0x0) return SPXERR_NULL_POINTER;
  struct wcserr **err = &(spx->err);

  int haverest = 1;
  if (restfrq == 0.0) {
    if (restwav == 0.0) {
      // No line rest frequency supplied.
      haverest = 0;

      // Temporarily set a dummy value for conversions.
      spx->restwav = 1.0;
    } else {
      spx->restwav = restwav;
    }
    spx->restfrq = C/spx->restwav;

  } else {
    spx->restfrq = restfrq;
    spx->restwav = C/restfrq;
  }

  spx->err = 0x0;

  // Convert to frequency.
  spx->wavetype = 0;
  spx->velotype = 0;
  if (strcmp(type, "FREQ") == 0) {
    if (spec == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable: frequency == 0");
    }
    spx->freq = spec;
    spx->wavetype = 1;

  } else if (strcmp(type, "AFRQ") == 0) {
    if (spec == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable: frequency == 0");
    }
    spx->freq = spec/(2.0*PI);
    spx->wavetype = 1;

  } else if (strcmp(type, "ENER") == 0) {
    if (spec == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable: frequency == 0");
    }
    spx->freq = spec/h;
    spx->wavetype = 1;

  } else if (strcmp(type, "WAVN") == 0) {
    if (spec == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable: frequency == 0");
    }
    spx->freq = spec*C;
    spx->wavetype = 1;

  } else if (strcmp(type, "VRAD") == 0) {
    spx->freq = spx->restfrq*(1.0 - spec/C);
    spx->velotype = 1;

  } else if (strcmp(type, "WAVE") == 0) {
    if (spec == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable: frequency == 0");
    }
    spx->freq = C/spec;
    spx->wavetype = 1;

  } else if (strcmp(type, "VOPT") == 0) {
    double s = 1.0 + spec/C;
    if (s == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable");
    }
    spx->freq = spx->restfrq/s;
    spx->velotype = 1;

  } else if (strcmp(type, "ZOPT") == 0) {
    double s = 1.0 + spec;
    if (s == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable");
    }
    spx->freq = spx->restfrq/s;
    spx->velotype = 1;

  } else if (strcmp(type, "AWAV") == 0) {
    if (spec == 0.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable");
    }

    double n, s;
    s = 1.0/spec;
    s *= s;
    n  =   2.554e8 / (0.41e14 - s);
    n += 294.981e8 / (1.46e14 - s);
    n += 1.000064328;
    spx->freq = C/(spec*n);
    spx->wavetype = 1;

  } else if (strcmp(type, "VELO") == 0) {
    double beta = spec/C;
    if (fabs(beta) == 1.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable");
    }
    spx->freq = spx->restfrq*(1.0 - beta)/sqrt(1.0 - beta*beta);
    spx->velotype = 1;

  } else if (strcmp(type, "BETA") == 0) {
    if (fabs(spec) == 1.0) {
      return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_VAR),
        "Invalid spectral variable");
    }
    spx->freq = spx->restfrq*(1.0 - spec)/sqrt(1.0 - spec*spec);
    spx->velotype = 1;

  } else {
    // Unrecognized type.
    return wcserr_set(WCSERR_SET(SPXERR_BAD_SPEC_PARAMS),
      "Unrecognized spectral type '%s'", type);
  }


  // Convert frequency to the other spectral types.
  double n, s, t, u;
  n = 1.0;
  for (int k = 0; k < 4; k++) {
    s = n*spx->freq/C;
    s *= s;
    t = 0.41e14 - s;
    u = 1.46e14 - s;
    n = 1.000064328 + (2.554e8/t + 294.981e8/u);
  }

  double dwaveawav = n - 2.0*s*(2.554e8/(t*t) + 294.981e8/(u*u));

  s = spx->freq/spx->restfrq;

  spx->ener = spx->freq*h;
  spx->afrq = spx->freq*(2.0*PI);
  spx->wavn = spx->freq/C;
  spx->vrad = C*(1.0 - s);
  spx->wave = C/spx->freq;
  spx->awav = spx->wave/n;
  spx->vopt = C*(1.0/s - 1.0);
  spx->zopt = spx->vopt/C;
  spx->velo = C*(1.0 - s*s)/(1.0 + s*s);
  spx->beta = spx->velo/C;

  // Compute the required derivatives.
  double gamma = 1.0/sqrt(1.0 - spx->beta*spx->beta);

  spx->dfreqafrq = 1.0/(2.0*PI);
  spx->dafrqfreq = 1.0/spx->dfreqafrq;

  spx->dfreqener = 1.0/h;
  spx->denerfreq = 1.0/spx->dfreqener;

  spx->dfreqwavn = C;
  spx->dwavnfreq = 1.0/spx->dfreqwavn;

  spx->dfreqvrad = -spx->restfrq/C;
  spx->dvradfreq = 1.0/spx->dfreqvrad;

  spx->dfreqwave = -spx->freq/spx->wave;
  spx->dwavefreq = 1.0/spx->dfreqwave;

  spx->dfreqawav = spx->dfreqwave * dwaveawav;
  spx->dawavfreq = 1.0/spx->dfreqawav;

  spx->dfreqvelo = -gamma*spx->restfrq/(C + spx->velo);
  spx->dvelofreq = 1.0/spx->dfreqvelo;

  spx->dwavevopt = spx->restwav/C;
  spx->dvoptwave = 1.0/spx->dwavevopt;

  spx->dwavezopt = spx->restwav;
  spx->dzoptwave = 1.0/spx->dwavezopt;

  spx->dwaveawav = dwaveawav;
  spx->dawavwave = 1.0/spx->dwaveawav;

  spx->dwavevelo = gamma*spx->restwav/(C - spx->velo);
  spx->dvelowave = 1.0/spx->dwavevelo;

  spx->dawavvelo = spx->dwavevelo/dwaveawav;
  spx->dveloawav = 1.0/spx->dawavvelo;

  spx->dvelobeta = C;
  spx->dbetavelo = 1.0/spx->dvelobeta;


  // Reset values if no line rest frequency was supplied.
  if (haverest) {
    spx->wavetype = 1;
    spx->velotype = 1;

  } else {
    spx->restfrq = 0.0;
    spx->restwav = 0.0;

    if (!spx->wavetype) {
      // Don't have wave characteristic types.
      spx->freq = 0.0;
      spx->afrq = 0.0;
      spx->ener = 0.0;
      spx->wavn = 0.0;
      spx->wave = 0.0;
      spx->awav = 0.0;

      spx->dfreqwave = 0.0;
      spx->dwavefreq = 0.0;

      spx->dfreqawav = 0.0;
      spx->dawavfreq = 0.0;

      spx->dwaveawav = 0.0;
      spx->dawavwave = 0.0;

    } else {
      // Don't have velocity types.
      spx->vrad = 0.0;
      spx->vopt = 0.0;
      spx->zopt = 0.0;
      spx->velo = 0.0;
      spx->beta = 0.0;
    }

    spx->dfreqvrad = 0.0;
    spx->dvradfreq = 0.0;

    spx->dfreqvelo = 0.0;
    spx->dvelofreq = 0.0;

    spx->dwavevopt = 0.0;
    spx->dvoptwave = 0.0;

    spx->dwavezopt = 0.0;
    spx->dzoptwave = 0.0;

    spx->dwavevelo = 0.0;
    spx->dvelowave = 0.0;

    spx->dawavvelo = 0.0;
    spx->dveloawav = 0.0;
  }

  return 0;
}

//----------------------------------------------------------------------------

int spxperr(const struct spxprm *spx, const char *prefix)

{
  if (spx == 0x0) return SPXERR_NULL_POINTER;

  if (spx->err) {
    wcserr_prt(spx->err, prefix);
  }

  return 0;
}


/*============================================================================
*   Conversions between frequency and vacuum wavelength.
*===========================================================================*/

int freqwave(
  double dummy,
  int nfreq,
  int sfreq,
  int swave,
  const double freq[],
  double wave[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *freqp = freq;
  double *wavep = wave;
  int *statp = stat;
  int status = 0;
  for (int ifreq = 0; ifreq < nfreq; ifreq++) {
    if (*freqp != 0.0) {
      *wavep = C/(*freqp);
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = SPXERR_BAD_INSPEC_COORD;
    }

    freqp += sfreq;
    wavep += swave;
  }

  return status;
}

//----------------------------------------------------------------------------

int wavefreq(
  double dummy,
  int nwave,
  int swave,
  int sfreq,
  const double wave[],
  double freq[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *wavep = wave;
  double *freqp = freq;
  int *statp = stat;
  int status = 0;
  for (int iwave = 0; iwave < nwave; iwave++) {
    if (*wavep != 0.0) {
      *freqp = C/(*wavep);
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = SPXERR_BAD_INSPEC_COORD;
    }

    wavep += swave;
    freqp += sfreq;
  }

  return status;
}

/*============================================================================
*   Conversions between frequency and air wavelength.
*===========================================================================*/

int freqawav(
  double dummy,
  int nfreq,
  int sfreq,
  int sawav,
  const double freq[],
  double awav[],
  int stat[])

{
  int status;
  if ((status = freqwave(dummy, nfreq, sfreq, sawav, freq, awav, stat))) {
    return status;
  }

  return waveawav(dummy, nfreq, sawav, sawav, awav, awav, stat);
}

//----------------------------------------------------------------------------

int awavfreq(
  double dummy,
  int nawav,
  int sawav,
  int sfreq,
  const double awav[],
  double freq[],
  int stat[])

{
  int status;
  if ((status = awavwave(dummy, nawav, sawav, sfreq, awav, freq, stat))) {
    return status;
  }

  return wavefreq(dummy, nawav, sfreq, sfreq, freq, freq, stat);
}

/*============================================================================
*   Conversions between frequency and relativistic velocity.
*===========================================================================*/

int freqvelo(
  double restfrq,
  int nfreq,
  int sfreq,
  int svelo,
  const double freq[],
  double velo[],
  int stat[])

{
  double r = restfrq*restfrq;

  const double *freqp = freq;
  double *velop = velo;
  int *statp = stat;
  for (int ifreq = 0; ifreq < nfreq; ifreq++) {
    double s = *freqp * *freqp;
    *velop = C*(r - s)/(r + s);
    *(statp++) = 0;

    freqp += sfreq;
    velop += svelo;
  }

  return 0;
}

//----------------------------------------------------------------------------

int velofreq(
  double restfrq,
  int nvelo,
  int svelo,
  int sfreq,
  const double velo[],
  double freq[],
  int stat[])

{
  const double *velop = velo;
  double *freqp = freq;
  int *statp = stat;
  int status = 0;
  for (int ivelo = 0; ivelo < nvelo; ivelo++) {
    double s = C + *velop;
    if (s != 0.0) {
      *freqp = restfrq*sqrt((C - *velop)/s);
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = SPXERR_BAD_INSPEC_COORD;
    }

    velop += svelo;
    freqp += sfreq;
  }

  return status;
}

/*============================================================================
*   Conversions between vacuum wavelength and air wavelength.
*===========================================================================*/

int waveawav(
  double dummy,
  int nwave,
  int swave,
  int sawav,
  const double wave[],
  double awav[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *wavep = wave;
  double *awavp = awav;
  int *statp = stat;
  int status = 0;
  for (int iwave = 0; iwave < nwave; iwave++) {
    if (*wavep != 0.0) {
      double n = 1.0;
      for (int k = 0; k < 4; k++) {
        double s  = n/(*wavep);
        s *= s;
        n  =   2.554e8 / (0.41e14 - s);
        n += 294.981e8 / (1.46e14 - s);
        n += 1.000064328;
      }

      *awavp = (*wavep)/n;
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = SPXERR_BAD_INSPEC_COORD;
    }

    wavep += swave;
    awavp += sawav;
  }

  return status;
}

//----------------------------------------------------------------------------

int awavwave(
  double dummy,
  int nawav,
  int sawav,
  int swave,
  const double awav[],
  double wave[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *awavp = awav;
  double *wavep = wave;
  int *statp = stat;
  int status = 0;
  for (int iawav = 0; iawav < nawav; iawav++) {
    if (*awavp != 0.0) {
      double n, s;
      s = 1.0/(*awavp);
      s *= s;
      n  =   2.554e8 / (0.41e14 - s);
      n += 294.981e8 / (1.46e14 - s);
      n += 1.000064328;
      *wavep = (*awavp)*n;
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = SPXERR_BAD_INSPEC_COORD;
    }

    awavp += sawav;
    wavep += swave;
  }

  return status;
}

/*============================================================================
*   Conversions between vacuum wavelength and relativistic velocity.
*===========================================================================*/

int wavevelo(
  double restwav,
  int nwave,
  int swave,
  int svelo,
  const double wave[],
  double velo[],
  int stat[])

{
  double r = restwav*restwav;

  const double *wavep = wave;
  double *velop = velo;
  int *statp = stat;
  for (int iwave = 0; iwave < nwave; iwave++) {
    double s = *wavep * *wavep;
    *velop = C*(s - r)/(s + r);
    *(statp++) = 0;

    wavep += swave;
    velop += svelo;
  }

  return 0;
}

//----------------------------------------------------------------------------

int velowave(
  double restwav,
  int nvelo,
  int svelo,
  int swave,
  const double velo[],
  double wave[],
  int stat[])

{
  const double *velop = velo;
  double *wavep = wave;
  int *statp = stat;
  int status = 0;
  for (int ivelo = 0; ivelo < nvelo; ivelo++) {
    double s = C - *velop;
    if (s != 0.0) {
      *wavep = restwav*sqrt((C + *velop)/s);
      *(statp++) = 0;
    } else {
      *(statp++) = 1;
      status = SPXERR_BAD_INSPEC_COORD;
    }

    velop += svelo;
    wavep += swave;
  }

  return status;
}

/*============================================================================
*   Conversions between air wavelength and relativistic velocity.
*===========================================================================*/

int awavvelo(
  double dummy,
  int nawav,
  int sawav,
  int svelo,
  const double awav[],
  double velo[],
  int stat[])

{
  int status;
  if ((status = awavwave(dummy, nawav, sawav, svelo, awav, velo, stat))) {
    return status;
  }

  return wavevelo(dummy, nawav, svelo, svelo, velo, velo, stat);
}

//----------------------------------------------------------------------------

int veloawav(
  double dummy,
  int nvelo,
  int svelo,
  int sawav,
  const double velo[],
  double awav[],
  int stat[])

{
  int status;
  if ((status = velowave(dummy, nvelo, svelo, sawav, velo, awav, stat))) {
    return status;
  }

  return waveawav(dummy, nvelo, sawav, sawav, awav, awav, stat);
}

/*============================================================================
*   Conversions between frequency and angular frequency.
*===========================================================================*/

int freqafrq(
  double dummy,
  int nfreq,
  int sfreq,
  int safrq,
  const double freq[],
  double afrq[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *freqp = freq;
  double *afrqp = afrq;
  int *statp = stat;
  for (int ifreq = 0; ifreq < nfreq; ifreq++) {
    *afrqp = (*freqp)*(2.0*PI);
    *(statp++) = 0;

    freqp += sfreq;
    afrqp += safrq;
  }

  return 0;
}

//----------------------------------------------------------------------------

int afrqfreq(
  double dummy,
  int nafrq,
  int safrq,
  int sfreq,
  const double afrq[],
  double freq[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *afrqp = afrq;
  double *freqp = freq;
  int *statp = stat;
  for (int iafrq = 0; iafrq < nafrq; iafrq++) {
    *freqp = (*afrqp)/(2.0*PI);
    *(statp++) = 0;

    afrqp += safrq;
    freqp += sfreq;
  }

  return 0;
}

/*============================================================================
*   Conversions between frequency and energy.
*===========================================================================*/

int freqener(
  double dummy,
  int nfreq,
  int sfreq,
  int sener,
  const double freq[],
  double ener[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *freqp = freq;
  double *enerp = ener;
  int *statp = stat;
  for (int ifreq = 0; ifreq < nfreq; ifreq++) {
    *enerp = (*freqp)*h;
    *(statp++) = 0;

    freqp += sfreq;
    enerp += sener;
  }

  return 0;
}

//----------------------------------------------------------------------------

int enerfreq(
  double dummy,
  int nener,
  int sener,
  int sfreq,
  const double ener[],
  double freq[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *enerp = ener;
  double *freqp = freq;
  int *statp = stat;
  for (int iener = 0; iener < nener; iener++) {
    *freqp = (*enerp)/h;
    *(statp++) = 0;

    enerp += sener;
    freqp += sfreq;
  }

  return 0;
}

/*============================================================================
*   Conversions between frequency and wave number.
*===========================================================================*/

int freqwavn(
  double dummy,
  int nfreq,
  int sfreq,
  int swavn,
  const double freq[],
  double wavn[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *freqp = freq;
  double *wavnp = wavn;
  int *statp = stat;
  for (int ifreq = 0; ifreq < nfreq; ifreq++) {
    *wavnp = (*freqp)/C;
    *(statp++) = 0;

    freqp += sfreq;
    wavnp += swavn;
  }

  return 0;
}

//----------------------------------------------------------------------------

int wavnfreq(
  double dummy,
  int nwavn,
  int swavn,
  int sfreq,
  const double wavn[],
  double freq[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *wavnp = wavn;
  double *freqp = freq;
  int *statp = stat;
  for (int iwavn = 0; iwavn < nwavn; iwavn++) {
    *freqp = (*wavnp)*C;
    *(statp++) = 0;

    wavnp += swavn;
    freqp += sfreq;
  }

  return 0;
}

/*============================================================================
*   Conversions between frequency and radio velocity.
*===========================================================================*/

int freqvrad(
  double restfrq,
  int nfreq,
  int sfreq,
  int svrad,
  const double freq[],
  double vrad[],
  int stat[])

{
  if (restfrq == 0.0) {
    return SPXERR_BAD_SPEC_PARAMS;
  }
  double r = C/restfrq;

  const double *freqp = freq;
  double *vradp = vrad;
  int *statp = stat;
  for (int ifreq = 0; ifreq < nfreq; ifreq++) {
    *vradp = r*(restfrq - *freqp);
    *(statp++) = 0;

    freqp += sfreq;
    vradp += svrad;
  }

  return 0;
}

//----------------------------------------------------------------------------

int vradfreq(
  double restfrq,
  int nvrad,
  int svrad,
  int sfreq,
  const double vrad[],
  double freq[],
  int stat[])

{
  double r = restfrq/C;

  const double *vradp = vrad;
  double *freqp = freq;
  int *statp = stat;
  for (int ivrad = 0; ivrad < nvrad; ivrad++) {
    *freqp = r*(C - *vradp);
    *(statp++) = 0;
    vradp += svrad;
    freqp += sfreq;
  }

  return 0;
}

/*============================================================================
*   Conversions between vacuum wavelength and optical velocity.
*===========================================================================*/

int wavevopt(
  double restwav,
  int nwave,
  int swave,
  int svopt,
  const double wave[],
  double vopt[],
  int stat[])

{
  if (restwav == 0.0) {
    return SPXERR_BAD_SPEC_PARAMS;
  }
  double r = C/restwav;

  const double *wavep = wave;
  double *voptp = vopt;
  int *statp = stat;
  for (int iwave = 0; iwave < nwave; iwave++) {
    *voptp = r*(*wavep) - C;
    *(statp++) = 0;
    wavep += swave;
    voptp += svopt;
  }

  return 0;
}

//----------------------------------------------------------------------------

int voptwave(
  double restwav,
  int nvopt,
  int svopt,
  int swave,
  const double vopt[],
  double wave[],
  int stat[])

{
  double r = restwav/C;

  const double *voptp = vopt;
  double *wavep = wave;
  int *statp = stat;
  for (int ivopt = 0; ivopt < nvopt; ivopt++) {
    *wavep = r*(C + *voptp);
    *(statp++) = 0;
    voptp += svopt;
    wavep += swave;
  }

  return 0;
}

/*============================================================================
*   Conversions between vacuum wavelength and redshift.
*===========================================================================*/

int wavezopt(
  double restwav,
  int nwave,
  int swave,
  int szopt,
  const double wave[],
  double zopt[],
  int stat[])

{
  if (restwav == 0.0) {
    return SPXERR_BAD_SPEC_PARAMS;
  }
  double r = 1.0/restwav;

  const double *wavep = wave;
  double *zoptp = zopt;
  int *statp = stat;
  for (int iwave = 0; iwave < nwave; iwave++) {
    *zoptp = r*(*wavep) - 1.0;
    *(statp++) = 0;
    wavep += swave;
    zoptp += szopt;
  }

  return 0;
}

//----------------------------------------------------------------------------

int zoptwave(
  double restwav,
  int nzopt,
  int szopt,
  int swave,
  const double zopt[],
  double wave[],
  int stat[])

{
  const double *zoptp = zopt;
  double *wavep = wave;
  int *statp = stat;
  for (int izopt = 0; izopt < nzopt; izopt++) {
    *wavep = restwav*(1.0 + *zoptp);
    *(statp++) = 0;
    zoptp += szopt;
    wavep += swave;
  }

  return 0;
}

/*============================================================================
*   Conversions between relativistic velocity and beta (= v/c).
*===========================================================================*/

int velobeta(
  double dummy,
  int nvelo,
  int svelo,
  int sbeta,
  const double velo[],
  double beta[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *velop = velo;
  double *betap = beta;
  int *statp = stat;
  for (int ivelo = 0; ivelo < nvelo; ivelo++) {
    *betap = (*velop)/C;
    *(statp++) = 0;

    velop += svelo;
    betap += sbeta;
  }

  return 0;
}

//----------------------------------------------------------------------------

int betavelo(
  double dummy,
  int nbeta,
  int sbeta,
  int svelo,
  const double beta[],
  double velo[],
  int stat[])

{
  // Avert nuisance compiler warnings about unused parameters.
  (void)dummy;

  const double *betap = beta;
  double *velop = velo;
  int *statp = stat;
  for (int ibeta = 0; ibeta < nbeta; ibeta++) {
    *velop = (*betap)*C;
    *(statp++) = 0;

    betap += sbeta;
    velop += svelo;
  }

  return 0;
}
