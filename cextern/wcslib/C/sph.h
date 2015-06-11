/*============================================================================

  WCSLIB 5.3 - an implementation of the FITS WCS standard.
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
  $Id: sph.h,v 5.3 2015/04/21 02:50:51 mcalabre Exp $
*=============================================================================
*
* WCSLIB 5.3 - C routines that implement the spherical coordinate
* transformations used by the FITS World Coordinate System (WCS) standard.
* Refer to
*
*   "Representations of world coordinates in FITS",
*   Greisen, E.W., & Calabretta, M.R. 2002, A&A, 395, 1061 (Paper I)
*
*   "Representations of celestial coordinates in FITS",
*   Calabretta, M.R., & Greisen, E.W. 2002, A&A, 395, 1077 (Paper II)
*
* Refer to the README file provided with WCSLIB for an overview of the
* library.
*
*
* Summary of the sph routines
* ---------------------------
* The WCS spherical coordinate transformations are implemented via separate
* functions, sphx2s() and sphs2x(), for the transformation in each direction.
*
* A utility function, sphdpa(), computes the angular distances and position
* angles from a given point on the sky to a number of other points.  sphpad()
* does the complementary operation - computes the coordinates of points offset
* by the given angular distances and position angles from a given point on the
* sky.
*
*
* sphx2s() - Rotation in the pixel-to-world direction
* ---------------------------------------------------
* sphx2s() transforms native coordinates of a projection to celestial
* coordinates.
*
* Given:
*   eul       const double[5]
*                       Euler angles for the transformation:
*                         0: Celestial longitude of the native pole [deg].
*                         1: Celestial colatitude of the native pole, or
*                            native colatitude of the celestial pole [deg].
*                         2: Native longitude of the celestial pole [deg].
*                         3: cos(eul[1])
*                         4: sin(eul[1])
*
*   nphi,
*   ntheta    int       Vector lengths.
*
*   spt,sxy   int       Vector strides.
*
*   phi,theta const double[]
*                       Longitude and latitude in the native coordinate
*                       system of the projection [deg].
*
* Returned:
*   lng,lat   double[]  Celestial longitude and latitude [deg].  These may
*                       refer to the same storage as phi and theta
*                       respectively.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
*
* sphs2x() - Rotation in the world-to-pixel direction
* ---------------------------------------------------
* sphs2x() transforms celestial coordinates to the native coordinates of a
* projection.
*
* Given:
*   eul       const double[5]
*                       Euler angles for the transformation:
*                         0: Celestial longitude of the native pole [deg].
*                         1: Celestial colatitude of the native pole, or
*                            native colatitude of the celestial pole [deg].
*                         2: Native longitude of the celestial pole [deg].
*                         3: cos(eul[1])
*                         4: sin(eul[1])
*
*   nlng,nlat int       Vector lengths.
*
*   sll,spt   int       Vector strides.
*
*   lng,lat   const double[]
*                       Celestial longitude and latitude [deg].
*
* Returned:
*   phi,theta double[]  Longitude and latitude in the native coordinate system
*                       of the projection [deg].  These may refer to the same
*                       storage as lng and lat respectively.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
*
* sphdpa() - Compute angular distance and position angle
* ------------------------------------------------------
* sphdpa() computes the angular distance and generalized position angle (see
* notes) from a "reference" point to a number of "field" points on the sphere.
* The points must be specified consistently in any spherical coordinate
* system.
*
* sphdpa() is complementary to sphpad().
*
* Given:
*   nfield    int       The number of field points.
*
*   lng0,lat0 double    Spherical coordinates of the reference point [deg].
*
*   lng,lat   const double[]
*                       Spherical coordinates of the field points [deg].
*
* Returned:
*   dist,pa   double[]  Angular distances and position angles [deg].  These
*                       may refer to the same storage as lng and lat
*                       respectively.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
* Notes:
*   sphdpa() uses sphs2x() to rotate coordinates so that the reference point
*   is at the north pole of the new system with the north pole of the old
*   system at zero longitude in the new.  The Euler angles required by
*   sphs2x() for this rotation are
*
=     eul[0] = lng0;
=     eul[1] = 90.0 - lat0;
=     eul[2] =  0.0;
*
*   The angular distance and generalized position angle are readily obtained
*   from the longitude and latitude of the field point in the new system.
*   This applies even if the reference point is at one of the poles, in which
*   case the "position angle" returned is as would be computed for a reference
*   point at (lng0,+90-epsilon) or (lng0,-90+epsilon), in the limit as epsilon
*   goes to zero.
*
*   It is evident that the coordinate system in which the two points are
*   expressed is irrelevant to the determination of the angular separation
*   between the points.  However, this is not true of the generalized position
*   angle.
*
*   The generalized position angle is here defined as the angle of
*   intersection of the great circle containing the reference and field points
*   with that containing the reference point and the pole.  It has its normal
*   meaning when the the reference and field points are specified in
*   equatorial coordinates (right ascension and declination).
*
*   Interchanging the reference and field points changes the position angle in
*   a non-intuitive way (because the sum of the angles of a spherical triangle
*   normally exceeds 180 degrees).
*
*   The position angle is undefined if the reference and field points are
*   coincident or antipodal.  This may be detected by checking for a distance
*   of 0 or 180 degrees (within rounding tolerance).  sphdpa() will return an
*   arbitrary position angle in such circumstances.
*
*
* sphpad() - Compute field points offset from a given point
* ---------------------------------------------------------
* sphpad() computes the coordinates of a set of points that are offset by the
* specified angular distances and position angles from a given "reference"
* point on the sky.  The distances and position angles must be specified
* consistently in any spherical coordinate system.
*
* sphpad() is complementary to sphdpa().
*
* Given:
*   nfield    int       The number of field points.
*
*   lng0,lat0 double    Spherical coordinates of the reference point [deg].
*
*   dist,pa   const double[]
*                       Angular distances and position angles [deg].
*
* Returned:
*   lng,lat   double[]  Spherical coordinates of the field points [deg].
*                       These may refer to the same storage as dist and pa
*                       respectively.
*
* Function return value:
*             int       Status return value:
*                         0: Success.
*
* Notes:
*   sphpad() is implemented analogously to sphdpa() although using sphx2s()
*   for the inverse transformation.  In particular, when the reference point
*   is at one of the poles, "position angle" is interpreted as though the
*   reference point was at (lng0,+90-epsilon) or (lng0,-90+epsilon), in the
*   limit as epsilon goes to zero.
*
*   Applying sphpad() with the distances and position angles computed by
*   sphdpa() should return the original field points.
*
*===========================================================================*/

#ifndef WCSLIB_SPH
#define WCSLIB_SPH

#ifdef __cplusplus
extern "C" {
#endif


int sphx2s(const double eul[5], int nphi, int ntheta, int spt, int sxy,
           const double phi[], const double theta[],
           double lng[], double lat[]);

int sphs2x(const double eul[5], int nlng, int nlat, int sll , int spt,
           const double lng[], const double lat[],
           double phi[], double theta[]);

int sphdpa(int nfield, double lng0, double lat0,
           const double lng[], const double lat[],
           double dist[], double pa[]);

int sphpad(int nfield, double lng0, double lat0,
           const double dist[], const double pa[],
           double lng[], double lat[]);


#ifdef __cplusplus
}
#endif

#endif /* WCSLIB_SPH */
