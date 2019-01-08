#include "erfa.h"

void eraPv2pav(double pv[2][3], double p[3], double v[3])
/*
**  - - - - - - - - -
**   e r a P v 2 p a v
**  - - - - - - - - -
**
**  Extend a p-vector to a pv-vector by appending a zero velocity.
**
**  Given:
**     pv       double[2][3]    pv-vector
**
**  Returned:
**     p        double[3]       p-vector
**     v        double[3]       v-vector
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013-2017, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraCp(pv[0], p);
   eraCp(pv[1], v);

   return;

}
