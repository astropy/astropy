:orphan:

==================================
FLRW Cosmology Method Cheat Sheet
==================================
 A class describing an isotropic and homogeneous (Friedmann-Lemaitre-Robertson-Walker) cosmology.

 Refer to `FLRW documentation for detailed information <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html/>`_.


Class Definition
+++++++++++++++++

``class astropy.cosmology.FLRW(H0, Om0, Ode0, Tcmb0=0, Neff=3.04, m_nu=<Quantity 0. eV>, Ob0=None, name=None)``


 =====================================================================================================================================================================    ==============================================================================================================================================
  Method                                                                                                                                                                   Summary
 =====================================================================================================================================================================    ==============================================================================================================================================
 `H(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.H>`_	                                                                     Hubble parameter (km/s/Mpc) at redshift z.
 `Ob(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Ob>`_	                                                                   Return the density parameter for baryonic matter at redshift z.
 `Ode(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Ode>`_	                                                                 Return the density parameter for dark energy at redshift z.
 `Odm(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Odm>`_	                                                                 Return the density parameter for dark matter at redshift z.
 `Ogamma(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Ogamma>`_	                                                           Return the density parameter for photons at redshift z.
 `Ok(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Ok>`_	                                                                   Return the equivalent density parameter for curvature at redshift z.
 `Om(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Om>`_	                                                                   Return the density parameter for non-relativistic matter at redshift z.
 `Onu(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Onu>`_	                                                                 Return the density parameter for neutrinos at redshift z.
 `Tcmb(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Tcmb>`_	                                                               Return the CMB temperature at redshift z.
 `Tnu(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.Tnu>`_	                                                                 Return the neutrino temperature at redshift z.
 `abs_distance_integrand(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.abs_distance_integrand>`_      	                     Integrand of the absorption distance.
 `absorption_distance(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.absorption_distance>`_	                                 Absorption distance at redshift z.
 `age(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.age>`_	                                                                 Age of the universe in Gyr at redshift z.
 `angular_diameter_distance(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.angular_diameter_distance>`_     	                 Angular diameter distance in Mpc at a given redshift.
 `angular_diameter_distance_z1z2(z1, z2) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.angular_diameter_distance_z1z2>`_ 	     Angular diameter distance between objects at 2 redshifts.
 `arcsec_per_kpc_comoving(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.arcsec_per_kpc_comoving>`_ 	                         Angular separation in arcsec corresponding to a comoving kpc at redshift z.
 `arcsec_per_kpc_proper(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.arcsec_per_kpc_proper>`_ 	                             Angular separation in arcsec corresponding to a proper kpc at redshift z.
 `clone(\*\*kwargs) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.clone>`_    	                                                 Returns a copy of this object, potentially with some changes.
 `comoving_distance(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.comoving_distance>`_ 	                                     Comoving line-of-sight distance in Mpc at a given redshift.
 `comoving_transverse_distance(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.comoving_transverse_distance>`_ 	               Comoving transverse distance in Mpc at a given redshift.
 `comoving_volume(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.comoving_volume>`_ 	                                         Comoving volume in cubic Mpc at redshift z.
 `critical_density(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.critical_density>`_ 	                                       Critical density in grams per cubic cm at redshift z.
 `de_density_scale(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.de_density_scale>`_ 	                                       Evaluates the redshift dependence of the dark energy density.
 `differential_comoving_volume(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.differential_comoving_volume>`_	               Differential comoving volume at redshift z.
 `distmod(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.distmod>`_	                                                         Distance modulus at redshift z.
 `efunc(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.efunc>`_	                                                             Function used to calculate H(z), the Hubble parameter.
 `inv_efunc(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.inv_efunc>`_	                                                     Inverse of efunc.
 `kpc_comoving_per_arcmin(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.kpc_comoving_per_arcmin>`_	                         Separation in transverse comoving kpc corresponding to an arcminute at redshift z.
 `kpc_proper_per_arcmin(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.kpc_proper_per_arcmin>`_	                             Separation in transverse proper kpc corresponding to an arcminute at redshift z.
 `lookback_distance(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.lookback_distance>`_ 	                                     The lookback distance is the light travel time distance to a given redshift.
 `lookback_time(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.lookback_time>`_ 	                                             Lookback time in Gyr to redshift z.
 `lookback_time_integrand(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.lookback_time_integrand>`_ 	                         Integrand of the lookback time.
 `luminosity_distance(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.luminosity_distance>`_ 	                                 Luminosity distance in Mpc at redshift z.
 `nu_relative_density(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.nu_relative_density>`_                                  Neutrino density function relative to the energy density in photons.
 `scale_factor(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.scale_factor>`_ 	                                               Scale factor at redshift z.
 `w(z) <http://docs.astropy.org/en/latest/api/astropy.cosmology.FLRW.html#astropy.cosmology.FLRW.w>`_ 	                                                                     The dark energy equation of state.
 =====================================================================================================================================================================    ==============================================================================================================================================
