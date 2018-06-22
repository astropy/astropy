#include "erfa.h"

void eraPav2pv(double p[3], double v[3], double pv[2][3])
/*
**  - - - - - - - - - -
**   e r a P a v 2 p v
**  - - - - - - - - - -
**
**  Extend a p-vector to a pv-vector by appending a zero velocity.
**
**  Given:
**     p        double[3]       p-vector
**     v        double[3]       v-vector
**
**  Returned:
**     pv       double[2][3]    pv-vector
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013-2017, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   eraCp(p, pv[0]);
   eraCp(v, pv[1]);

   return;

}
