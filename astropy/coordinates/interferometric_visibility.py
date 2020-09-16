"""
Provide the interferometric visibility frame, also known as the UVW frame, of radio astronomy to a precision
suitable for VLBI. The accuracy is comparable to CALC11.

In order to compute the coordinates in this frame we follow the approach of computing the proper time for a given
geodesic (which forms the w coordinate) and then computing the u and v coordinates from the directional derivatives of
proper time w.r.t the l and m coordinates.

w = proper-time
u = dw/dl
v = dw/dm

In order to perform VLBI-precise calcuations we must take into account the velocity of the local reference frame of each
geodesic origin.

TODO(@tikk3r, @joshuaalbert): derive the exact solution
TODO(@tikk3r, @joshuaalbert): Investigate setting up the frame using astropy.coordinates.SkyOffset
TODO(@tikk3r, @joshuaalbert): If possible, set up the frame using astropy.coordinates.SkyOffset
TODO(@tikk3r, @joshuaalbert): Set up CALC11 server to do testing
TODO(@tikk3r, @joshuaalbert): Investigate setting up CALC11 server as part of continuous integration
"""
pass