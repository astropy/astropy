/*============================================================================

  WCSLIB 5.2 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2015, Mark Calabretta

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

  Direct correspondence concerning WCSLIB to mark@calabretta.id.au

  Author: Mark Calabretta, Australia Telescope National Facility, CSIRO.
  http://www.atnf.csiro.au/people/Mark.Calabretta
  $Id: spx.h,v 5.2 2015/04/15 12:35:07 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.2 - C routines that implement the spectral coordinate systems
* recognized by the FITS World Coordinate System (WCS) standard.  Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (Paper I)
*
*   "Representations of spectral coordinates in FITS",
*   Greisen, E.W., Calabretta, M.R., Valdes, F.G., & Allen, S.L.
*   2006, A&A, 446, 747 (Paper III)
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the spx routines
* ---------------------------
* specx() is a scalar routine that, given one spectral variable (e.g.
* frequency), computes all the others (e.g. wavelength, velocity, etc.) plus
* the required derivatives of each with respect to the others.  The results
* are returned in the spxprm struct.
*
* The remaining routines are all vector conversions from one spectral
* variable to another.  The API of these functions only differ in whether the
* rest frequency or wavelength need be supplied.
*
* Non-linear:
*   - freqwave()    frequency              ->  vacuum wavelength
*   - wavefreq()    vacuum wavelength      ->  frequency
*
*   - freqawav()    frequency              ->  air wavelength
*   - awavfreq()    air wavelength         ->  frequency
*
*   - freqvelo()    frequency              ->  relativistic velocity
*   - velofreq()    relativistic velocity  ->  frequency
*
*   - waveawav()    vacuum wavelength      ->  air wavelength
*   - awavwave()    air wavelength         ->  vacuum wavelength
*
*   - wavevelo()    vacuum wavelength      ->  relativistic velocity
*   - velowave()    relativistic velocity  ->  vacuum wavelength
*
*   - awavvelo()    air wavelength         ->  relativistic velocity
*   - veloawav()    relativistic velocity  ->  air wavelength
*
* Linear:
*   - freqafrq()    frequency              ->  angular frequency
*   - afrqfreq()    angular frequency      ->  frequency
*
*   - freqener()    frequency              ->  energy
*   - enerfreq()    energy                 ->  frequency
*
*   - freqwavn()    frequency              ->  wave number
*   - wavnfreq()    wave number            ->  frequency
*
*   - freqvrad()    frequency              ->  radio velocity
*   - vradfreq()    radio velocity         ->  frequency
*
*   - wavevopt()    vacuum wavelength      ->  optical velocity
*   - voptwave()    optical velocity       ->  vacuum wavelength
*
*   - wavezopt()    vacuum wavelength      ->  redshift
*   - zoptwave()    redshift               ->  vacuum wavelength
*
*   - velobeta()    relativistic velocity  ->  beta (= v/c)
*   - betavelo()    beta (= v/c)           ->  relativistic velocity
*
* These are the workhorse routines, to be used for fast transformations.
* Conversions may be done "in place" by calling the routine with the output
* vector set to the input.
*
* Argument checking:
* ------------------
* The input spectral values are only checked for values that would result
* in floating point exceptions.  In particular, negative frequencies and
* wavelengths are allowed, as are velocities greater than the speed of
* light.  The same is true for the spectral parameters - rest frequency and
* wavelength.
*
* Accuracy:
* ---------
* No warranty is given for the accuracy of these routines (refer to the
* copyright notice); intending users must satisfy for themselves their
* adequacy for the intended purpose.  However, closure effectively to within
* double precision rounding error was demonstrated by test routine tspec.c
* which accompanies this software.
*
*
* specx() - Spectral cross conversions (scalar)
* ---------------------------------------------
* Given one spectral variable specx() computes all the others, plus the
* required derivatives of each with respect to the others.
*
* Given:
*   type      const char*
*                       The type of spectral variable given by spec, FREQ,
*                       AFRQ, ENER, WAVN, VRAD, WAVE, VOPT, ZOPT, AWAV, VELO,
*                       or BETA (case sensitive).
*
*   spec      double    The spectral variable given, in SI units.
*
*   restfrq,
*   restwav   double    Rest frequency [Hz] or rest wavelength in vacuo [m],
*                       only one of which need be given.  The other should be
*                       set to zero.  If both are zero, only a subset of the
*                       spectral variables can be computed, the remainder are
*                       set to zero.  Specifically, given one of FREQ, AFRQ,
*                       ENER, WAVN, WAVE, or AWAV the others can be computed
*                       without knowledge of the rest frequency.  Likewise,
*                       VRAD, VOPT, ZOPT, VELO, and BETA.
*
* Given and returned:
*   specs     struct spxprm*
*                       Data structure containing all spectral variables and
*                       their derivatives, in SI units.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         1: Null spxprm pointer passed.
*                         2: Invalid spectral parameters.
*                         3: Invalid spectral variable.
*
*                       For returns > 1, a detailed error message is set in
*                       spxprm::err if enabled, see wcserr_enable().
*
* freqafrq(), afrqfreq(), freqener(), enerfreq(), freqwavn(), wavnfreq(),
* freqwave(), wavefreq(), freqawav(), awavfreq(), waveawav(), awavwave(),
* velobeta(), and betavelo() implement vector conversions between wave-like
* or velocity-like spectral types (i.e. conversions that do not need the rest
* frequency or wavelength).  They all have the same API.
*
*
* freqafrq() - Convert frequency to angular frequency (vector)
* ------------------------------------------------------------
* freqafrq() converts frequency to angular frequency.
*
* Given:
*   param     double    Ignored.
*
*   nspec     int       Vector length.
*
*   instep,
*   outstep   int       Vector strides.
*
*   inspec    const double[]
*                       Input spectral variables, in SI units.
*
* Returned:
*   outspec   double[]  Output spectral variables, in SI units.
*
*   stat      int[]     Status return value for each vector element:
*                         0: Success.
*                         1: Invalid value of inspec.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters.
*                         4: One or more of the inspec coordinates were
*                            invalid, as indicated by the stat vector.
*
*
* freqvelo(), velofreq(), freqvrad(), and vradfreq() implement vector
* conversions between frequency and velocity spectral types.  They all have
* the same API.
*
*
* freqvelo() - Convert frequency to relativistic velocity (vector)
* ----------------------------------------------------------------
* freqvelo() converts frequency to relativistic velocity.
*
* Given:
*   param     double    Rest frequency [Hz].
*
*   nspec     int       Vector length.
*
*   instep,
*   outstep   int       Vector strides.
*
*   inspec    const double[]
*                       Input spectral variables, in SI units.
*
* Returned:
*   outspec   double[]  Output spectral variables, in SI units.
*
*   stat      int[]     Status return value for each vector element:
*                         0: Success.
*                         1: Invalid value of inspec.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters.
*                         4: One or more of the inspec coordinates were
*                            invalid, as indicated by the stat vector.
*
*
* wavevelo(), velowave(), awavvelo(), veloawav(), wavevopt(), voptwave(),
* wavezopt(), and zoptwave() implement vector conversions between wavelength
* and velocity spectral types.  They all have the same API.
*
*
* wavevelo() - Conversions between wavelength and velocity types (vector)
* -----------------------------------------------------------------------
* wavevelo() converts vacuum wavelength to relativistic velocity.
*
* Given:
*   param     double    Rest wavelength in vacuo [m].
*
*   nspec     int       Vector length.
*
*   instep,
*   outstep   int       Vector strides.
*
*   inspec    const double[]
*                       Input spectral variables, in SI units.
*
* Returned:
*   outspec   double[]  Output spectral variables, in SI units.
*
*   stat      int[]     Status return value for each vector element:
*                         0: Success.
*                         1: Invalid value of inspec.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*                         2: Invalid spectral parameters.
*                         4: One or more of the inspec coordinates were
*                            invalid, as indicated by the stat vector.
*
*
* spxprm struct - Spectral variables and their derivatives
* --------------------------------------------------------
* The spxprm struct contains the value of all spectral variables and their
* derivatives.   It is used solely by specx() which constructs it from
* information provided via its function arguments.
*
* This struct should be considered read-only, no members need ever be set nor
* should ever be modified by the user.
*
*   double restfrq
*     (Returned) Rest frequency [Hz].
*
*   double restwav
*     (Returned) Rest wavelength [m].
*
*   int wavetype
*     (Returned) True if wave types have been computed, and ...
*
*   int velotype
*     (Returned) ... true if velocity types have been computed; types are
*     defined below.
*
*     If one or other of spxprm::restfrq and spxprm::restwav is given
*     (non-zero) then all spectral variables may be computed.  If both are
*     given, restfrq is used.  If restfrq and restwav are both zero, only wave
*     characteristic xor velocity type spectral variables may be computed
*     depending on the variable given.   These flags indicate what is
*     available.
*
*   double freq
*     (Returned) Frequency [Hz] (wavetype).
*
*   double afrq
*     (Returned) Angular frequency [rad/s] (wavetype).
*
*   double ener
*     (Returned) Photon energy [J] (wavetype).
*
*   double wavn
*     (Returned) Wave number [/m] (wavetype).
*
*   double vrad
*     (Returned) Radio velocity [m/s] (velotype).
*
*   double wave
*     (Returned) Vacuum wavelength [m] (wavetype).
*
*   double vopt
*     (Returned) Optical velocity [m/s] (velotype).
*
*   double zopt
*     (Returned) Redshift [dimensionless] (velotype).
*
*   double awav
*     (Returned) Air wavelength [m] (wavetype).
*
*   double velo
*     (Returned) Relativistic velocity [m/s] (velotype).
*
*   double beta
*     (Returned) Relativistic beta [dimensionless] (velotype).
*
*   double dfreqafrq
*     (Returned) Derivative of frequency with respect to angular frequency
*     [/rad] (constant, = 1 / 2*pi), and ...
*   double dafrqfreq
*     (Returned) ... vice versa [rad] (constant, = 2*pi, always available).
*
*   double dfreqener
*     (Returned) Derivative of frequency with respect to photon energy
*     [/J/s] (constant, = 1/h), and ...
*   double denerfreq
*     (Returned) ... vice versa [Js] (constant, = h, Planck's constant,
*     always available).
*
*   double dfreqwavn
*     (Returned) Derivative of frequency with respect to wave number [m/s]
*     (constant, = c, the speed of light in vacuo), and ...
*   double dwavnfreq
*     (Returned) ... vice versa [s/m] (constant, = 1/c, always available).
*
*   double dfreqvrad
*     (Returned) Derivative of frequency with respect to radio velocity [/m],
*     and ...
*   double dvradfreq
*     (Returned) ... vice versa [m] (wavetype && velotype).
*
*   double dfreqwave
*     (Returned) Derivative of frequency with respect to vacuum wavelength
*     [/m/s], and ...
*   double dwavefreq
*     (Returned) ... vice versa [m s] (wavetype).
*
*   double dfreqawav
*     (Returned) Derivative of frequency with respect to air wavelength,
*     [/m/s], and ...
*   double dawavfreq
*     (Returned) ... vice versa [m s] (wavetype).
*
*   double dfreqvelo
*     (Returned) Derivative of frequency with respect to relativistic
*     velocity [/m], and ...
*   double dvelofreq
*     (Returned) ... vice versa [m] (wavetype && velotype).
*
*   double dwavevopt
*     (Returned) Derivative of vacuum wavelength with respect to optical
*     velocity [s], and ...
*   double dvoptwave
*     (Returned) ... vice versa [/s] (wavetype && velotype).
*
*   double dwavezopt
*     (Returned) Derivative of vacuum wavelength with respect to redshift [m],
*     and ...
*   double dzoptwave
*     (Returned) ... vice versa [/m] (wavetype && velotype).
*
*   double dwaveawav
*     (Returned) Derivative of vacuum wavelength with respect to air
*     wavelength [dimensionless], and ...
*   double dawavwave
*     (Returned) ... vice versa [dimensionless] (wavetype).
*
*   double dwavevelo
*     (Returned) Derivative of vacuum wavelength with respect to relativistic
*     velocity [s], and ...
*   double dvelowave
*     (Returned) ... vice versa [/s] (wavetype && velotype).
*
*   double dawavvelo
*     (Returned) Derivative of air wavelength with respect to relativistic
*     velocity [s], and ...
*   double dveloawav
*     (Returned) ... vice versa [/s] (wavetype && velotype).
*
*   double dvelobeta
*     (Returned) Derivative of relativistic velocity with respect to
*     relativistic beta [m/s] (constant, = c, the speed of light in vacuo),
*     and ...
*   double dbetavelo
*     (Returned) ... vice versa [s/m] (constant, = 1/c, always available).
*
*   struct wcserr *err
*     (Returned) If enabled, when an error status is returned, this struct
*     contains detailed information about the error, see wcserr_enable().
*
*   void *padding
*     (An unused variable inserted for alignment purposes only.)
*
* Global variable: const char *spx_errmsg[] - Status return messages
* ------------------------------------------------------------------
* Error messages to match the status value returned from each function.
*
*===========================================================================*/

#ifndef WCSLIB_SPEC
#define WCSLIB_SPEC

#ifdef __cplusplus
extern "C" {
#endif

extern const char *spx_errmsg[];

enum spx_errmsg {
  SPXERR_SUCCESS          = 0,	/* Success. */
  SPXERR_NULL_POINTER     = 1,	/* Null spxprm pointer passed. */
  SPXERR_BAD_SPEC_PARAMS  = 2,	/* Invalid spectral parameters. */
  SPXERR_BAD_SPEC_VAR     = 3,	/* Invalid spectral variable. */
  SPXERR_BAD_INSPEC_COORD = 4 	/* One or more of the inspec coordinates were
				   invalid. */
};

struct spxprm {
  double restfrq, restwav;	/* Rest frequency [Hz] and wavelength [m].  */

  int wavetype, velotype;	/* True if wave/velocity types have been    */
				/* computed; types are defined below.       */

  /* Spectral variables computed by specx().                                */
  /*------------------------------------------------------------------------*/
  double freq,			/* wavetype: Frequency [Hz].                */
         afrq,			/* wavetype: Angular frequency [rad/s].     */
         ener,			/* wavetype: Photon energy [J].             */
         wavn,			/* wavetype: Wave number [/m].              */
         vrad,			/* velotype: Radio velocity [m/s].          */
         wave,			/* wavetype: Vacuum wavelength [m].         */
         vopt,			/* velotype: Optical velocity [m/s].        */
         zopt,			/* velotype: Redshift.                      */
         awav,			/* wavetype: Air wavelength [m].            */
         velo,			/* velotype: Relativistic velocity [m/s].   */
         beta;			/* velotype: Relativistic beta.             */

  /* Derivatives of spectral variables computed by specx().                 */
  /*------------------------------------------------------------------------*/
  double dfreqafrq, dafrqfreq,	/* Constant, always available.              */
         dfreqener, denerfreq,	/* Constant, always available.              */
         dfreqwavn, dwavnfreq,	/* Constant, always available.              */
         dfreqvrad, dvradfreq,	/* wavetype && velotype.                    */
         dfreqwave, dwavefreq,	/* wavetype.                                */
         dfreqawav, dawavfreq,	/* wavetype.                                */
         dfreqvelo, dvelofreq,	/* wavetype && velotype.                    */
         dwavevopt, dvoptwave,	/* wavetype && velotype.                    */
         dwavezopt, dzoptwave,	/* wavetype && velotype.                    */
         dwaveawav, dawavwave,	/* wavetype.                                */
         dwavevelo, dvelowave,	/* wavetype && velotype.                    */
         dawavvelo, dveloawav,	/* wavetype && velotype.                    */
         dvelobeta, dbetavelo;	/* Constant, always available.              */

  /* Error handling                                                         */
  /*------------------------------------------------------------------------*/
  struct wcserr *err;

  /* Private                                                                */
  /*------------------------------------------------------------------------*/
  void   *padding;		/* (Dummy inserted for alignment purposes.) */
};

/* Size of the spxprm struct in int units, used by the Fortran wrappers. */
#define SPXLEN (sizeof(struct spxprm)/sizeof(int))


int specx(const char *type, double spec, double restfrq, double restwav,
          struct spxprm *specs);


/* For use in declaring function prototypes, e.g. in spcprm. */
#define SPX_ARGS double param, int nspec, int instep, int outstep, \
    const double inspec[], double outspec[], int stat[]

int freqafrq(SPX_ARGS);
int afrqfreq(SPX_ARGS);

int freqener(SPX_ARGS);
int enerfreq(SPX_ARGS);

int freqwavn(SPX_ARGS);
int wavnfreq(SPX_ARGS);

int freqwave(SPX_ARGS);
int wavefreq(SPX_ARGS);

int freqawav(SPX_ARGS);
int awavfreq(SPX_ARGS);

int waveawav(SPX_ARGS);
int awavwave(SPX_ARGS);

int velobeta(SPX_ARGS);
int betavelo(SPX_ARGS);


int freqvelo(SPX_ARGS);
int velofreq(SPX_ARGS);

int freqvrad(SPX_ARGS);
int vradfreq(SPX_ARGS);


int wavevelo(SPX_ARGS);
int velowave(SPX_ARGS);

int awavvelo(SPX_ARGS);
int veloawav(SPX_ARGS);

int wavevopt(SPX_ARGS);
int voptwave(SPX_ARGS);

int wavezopt(SPX_ARGS);
int zoptwave(SPX_ARGS);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_SPEC */
