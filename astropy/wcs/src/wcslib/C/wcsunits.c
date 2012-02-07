/*============================================================================

  WCSLIB 4.10 - an implementation of the FITS WCS standard.
  Copyright (C) 1995-2012, Mark Calabretta

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
  along with WCSLIB.  If not, see <http://www.gnu.org/licenses/>.

  Correspondence concerning WCSLIB may be directed to:
    Internet email: mcalabre@atnf.csiro.au
    Postal address: Dr. Mark Calabretta
                    Australia Telescope National Facility, CSIRO
                    PO Box 76
                    Epping NSW 1710
                    AUSTRALIA

  Author: Mark Calabretta, Australia Telescope National Facility
  http://www.atnf.csiro.au/~mcalabre/index.html
  $Id: wcsunits.c,v 4.10 2012/02/05 23:41:44 cal103 Exp $
*===========================================================================*/

#include <math.h>

#include "wcsunits.h"

/* Map status return value to message. */
const char *wcsunits_errmsg[] = {
  "Success",
  "Invalid numeric multiplier",
  "Dangling binary operator",
  "Invalid symbol in INITIAL context",
  "Function in invalid context",
  "Invalid symbol in EXPON context",
  "Unbalanced bracket",
  "Unbalanced parenthesis",
  "Consecutive binary operators",
  "Internal parser error",
  "Non-conformant unit specifications",
  "Non-conformant functions",
  "Potentially unsafe translation"};


/* Unit types. */
const char *wcsunits_types[] = {
  "plane angle",
  "solid angle",
  "charge",
  "mole",
  "temperature",
  "luminous intensity",
  "mass",
  "length",
  "time",
  "beam",
  "bin",
  "bit",
  "count",
  "stellar magnitude",
  "pixel",
  "solar ratio",
  "voxel"};

const char *wcsunits_units[] = {
  "degree",
  "steradian",
  "Coulomb",
  "mole",
  "Kelvin",
  "candela",
  "kilogram",
  "metre",
  "second",
  "", "", "", "", "", "", "", ""};

const char *wcsunits_funcs[] = {
  "none",
  "log",
  "ln",
  "exp"};

/*--------------------------------------------------------------------------*/

int wcsunits(
  const char have[],
  const char want[],
  double *scale,
  double *offset,
  double *power)

{
  return wcsunitse(
    have, want, scale, offset, power, 0x0);
}

/* : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : :  */

int wcsunitse(
  const char have[],
  const char want[],
  double *scale,
  double *offset,
  double *power,
  struct wcserr **err)

{
  static const char *function = "wcsunitse";

  int    func1, func2, i, status;
  double scale1, scale2, units1[WCSUNITS_NTYPE], units2[WCSUNITS_NTYPE];

  if ((status = wcsulexe(have, &func1, &scale1, units1, err))) {
    return status;
  }

  if ((status = wcsulexe(want, &func2, &scale2, units2, err))) {
    return status;
  }

  /* Check conformance. */
  for (i = 0; i < WCSUNITS_NTYPE; i++) {
    if (units1[i] != units2[i]) {
      return wcserr_set(WCSERR_SET(UNITSERR_BAD_UNIT_SPEC),
        "Mismatched units type '%s': have '%s', want '%s'",
        wcsunits_types[i], have, want);
    }
  }

  *scale  = 0.0;
  *offset = 0.0;
  *power  = 1.0;

  switch (func1) {
  case 0:
    /* No function. */
    if (func2) {
      return wcserr_set(WCSERR_SET(UNITSERR_BAD_FUNCS),
        "Mismatched unit functions: have '%s' (%s), want '%s' (%s)",
        have, wcsunits_funcs[func1], want, wcsunits_funcs[func2]);
    }

    *scale = scale1 / scale2;
    break;

  case 1:
    /* log(). */
    if (func2 == 1) {
      /* log(). */
      *scale  = 1.0;
      *offset = log10(scale1 / scale2);

    } else if (func2 == 2) {
      /* ln(). */
      *scale  = log(10.0);
      *offset = log(scale1 / scale2);

    } else {
      return wcserr_set(WCSERR_SET(UNITSERR_BAD_FUNCS),
        "Mismatched unit functions: have '%s' (%s), want '%s' (%s)",
        have, wcsunits_funcs[func1], want, wcsunits_funcs[func2]);
    }

    break;

  case 2:
    /* ln(). */
    if (func2 == 1) {
      /* log(). */
      *scale  = 1.0 / log(10.0);
      *offset = log(scale1 / scale2);

    } else if (func2 == 2) {
      /* ln(). */
      *scale  = 1.0;
      *offset = log(scale1 / scale2);

    } else {
      return wcserr_set(WCSERR_SET(UNITSERR_BAD_FUNCS),
        "Mismatched unit functions: have '%s' (%s), want '%s' (%s)",
        have, wcsunits_funcs[func1], want, wcsunits_funcs[func2]);
    }

    break;

  case 3:
    /* exp(). */
    if (func2 != 3) {
      return wcserr_set(WCSERR_SET(UNITSERR_BAD_FUNCS),
        "Mismatched unit functions: have '%s' (%s), want '%s' (%s)",
        have, wcsunits_funcs[func1], want, wcsunits_funcs[func2]);
    }

    *scale = 1.0;
    *power = scale1 / scale2;
    break;

  default:
    /* Internal parser error. */
    return wcserr_set(WCSERR_SET(UNITSERR_PARSER_ERROR),
      "Internal units parser error");
  }

  return 0;
}

/*--------------------------------------------------------------------------*/

int wcsutrn(int ctrl, char unitstr[])

{
  return wcsutrne(ctrl, unitstr, 0x0);
}

/*--------------------------------------------------------------------------*/

int wcsulex(const char unitstr[], int *func, double *scale, double units[])

{
  return wcsulexe(unitstr, func, scale, units, 0x0);
}
