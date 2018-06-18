#ifndef ERFAADDITIONSDEF
#define ERFAADDITIONSDEF

/*
**  - - - - - - - - - - - - - - - - -
**   e r f a _ a d d i t i o n s . h
**  - - - - - - - - - - - - - - - - -
**
**  A few extra routines which are particularly handy for constructing
**  pv vectors inside the coordinate transforms.
**
**  MHvK proposed these to Catherine Hohenkerk for inclusion in SOFA
**  on 2018-05-24, with the response suggesting this was reasonable and
**  might thus be done.
*/

/* Extra/PVMergeExtract */
void eraPav2pv(double p[3], double v[3], double pv[2][3]);
void eraPv2pav(double pv[2][3], double p[3], double v[3]);

#endif
