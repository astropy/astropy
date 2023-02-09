#cython: language_level=3, boundscheck=False
""" Cython inverse efuncs for cosmology integrals"""

cimport cython
from libc.math cimport exp, pow

## Inverse efunc methods for various dark energy subclasses
## These take only scalar arguments since that is what the integral
## routines give them.

## Implementation notes:
##  * Using a python list for nu_y seems to be faster than a ndarray,
##     given that nu_y generally has a small number of elements,
##     even when you turn off bounds checking, etc.
##  * Using pow(x, -0.5) is slightly faster than x**(-0.5) and
##    even more so than 1.0 / sqrt(x)
##  * Hardwiring in the p, 1/p, k, prefac values in nufunc is
##       nontrivially faster than declaring them with cdef

######### LambdaCDM
# No relativistic species
def lcdm_inv_efunc_norel(double z, double Om0, double Ode0, double Ok0):
  cdef double opz = 1.0 + z
  return pow(opz**2 * (opz * Om0 + Ok0) + Ode0, -0.5)

# Massless neutrinos
def lcdm_inv_efunc_nomnu(double z, double Om0, double Ode0, double Ok0,
    double Or0):
  cdef double opz = 1.0 + z
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 + Ode0, -0.5)

# With massive neutrinos
def lcdm_inv_efunc(double z, double Om0, double Ode0, double Ok0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 + Ode0, -0.5)

######## FlatLambdaCDM
# No relativistic species
def flcdm_inv_efunc_norel(double z, double Om0, double Ode0):
  return pow((1. + z)**3 * Om0 + Ode0, -0.5)

# Massless neutrinos
def flcdm_inv_efunc_nomnu(double z, double Om0, double Ode0, double Or0):
  cdef double opz = 1.0 + z
  return pow(opz**3 * (opz * Or0 + Om0) + Ode0, -0.5)

# With massive neutrinos
def flcdm_inv_efunc(double z, double Om0, double Ode0, double Ogamma0,
    double NeffPerNu, int nmasslessnu, list nu_y):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  return pow(opz**3 * (opz * Or0 + Om0) + Ode0, -0.5)

######## wCDM
# No relativistic species
def wcdm_inv_efunc_norel(double z, double Om0, double Ode0,
    double Ok0, double w0):
  cdef double opz = 1.0 + z
  return pow(opz**2 * (opz * Om0 + Ok0) +
            Ode0 * opz**(3. * (1.0 + w0)), -0.5)

# Massless neutrinos
def wcdm_inv_efunc_nomnu(double z, double Om0, double Ode0, double Ok0,
    double Or0, double w0):
  cdef double opz = 1.0 + z
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 +
          Ode0 * opz**(3. * (1.0 + w0)), -0.5)

# With massive neutrinos
def wcdm_inv_efunc(double z, double Om0, double Ode0, double Ok0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y, double w0):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 +
          Ode0 * opz**(3. * (1.0 + w0)), -0.5)

######## Flat wCDM
# No relativistic species
def fwcdm_inv_efunc_norel(double z, double Om0, double Ode0, double w0):
  cdef double opz = 1.0 + z
  return pow(opz**3 * Om0 + Ode0 * opz**(3. * (1.0 + w0)), -0.5)

# Massless neutrinos
def fwcdm_inv_efunc_nomnu(double z, double Om0, double Ode0,
    double Or0, double w0):
  cdef double opz = 1.0 + z
  return pow(opz**3 * (opz * Or0 + Om0) +
            Ode0 * opz**(3. * (1.0 + w0)), -0.5)

# With massive neutrinos
def fwcdm_inv_efunc(double z, double Om0, double Ode0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y, double w0):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  return pow(opz**3 * (opz * Or0 + Om0) + Ode0 * opz**(3. * (1.0 + w0)), -0.5)

######## w0waCDM
# No relativistic species
def w0wacdm_inv_efunc_norel(double z, double Om0, double Ode0, double Ok0,
    double w0, double wa):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1 + w0 + wa)) * exp(-3.0 * wa * z / opz)
  return pow(opz**2 * (opz * Om0 + Ok0) + Ode0 * Odescl, -0.5)

# Massless neutrinos
def w0wacdm_inv_efunc_nomnu(double z, double Om0, double Ode0, double Ok0,
    double Or0, double w0, double wa):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1 + w0 + wa)) * exp(-3.0 * wa * z / opz)
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 +
          Ode0 * Odescl, -0.5)

def w0wacdm_inv_efunc(double z, double Om0, double Ode0, double Ok0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y, double w0,
    double wa):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  cdef Odescl = opz**(3. * (1 + w0 + wa)) * exp(-3.0 * wa * z / opz)
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 + Ode0 * Odescl, -0.5)

######## Flatw0waCDM
# No relativistic species
def fw0wacdm_inv_efunc_norel(double z, double Om0, double Ode0,
    double w0, double wa):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1 + w0 + wa)) * exp(-3.0 * wa * z / opz)
  return pow(opz**3 * Om0 + Ode0 * Odescl, -0.5)

# Massless neutrinos
def fw0wacdm_inv_efunc_nomnu(double z, double Om0, double Ode0,
    double Or0, double w0, double wa):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1 + w0 + wa)) * exp(-3.0 * wa * z / opz)
  return pow((opz * Or0 + Om0) * opz**3 + Ode0 * Odescl, -0.5)

# With massive neutrinos
def fw0wacdm_inv_efunc(double z, double Om0, double Ode0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y, double w0,
    double wa):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  cdef Odescl = opz**(3. * (1 + w0 + wa)) * exp(-3.0 * wa * z / opz)
  return pow((opz * Or0 + Om0) * opz**3 + Ode0 * Odescl, -0.5)

######## wpwaCDM
# No relativistic species
def wpwacdm_inv_efunc_norel(double z, double Om0, double Ode0, double Ok0,
    double wp, double apiv, double wa):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1. + wp + apiv * wa)) * exp(-3. * wa * z / opz)
  return pow(opz**2 * (opz * Om0 + Ok0) + Ode0 * Odescl, -0.5)

# Massless neutrinos
def wpwacdm_inv_efunc_nomnu(double z, double Om0, double Ode0, double Ok0,
    double Or0, double wp, double apiv, double wa):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1. + wp + apiv * wa)) * exp(-3. * wa * z / opz)
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 +
          Ode0 * Odescl, -0.5)

# With massive neutrinos
def wpwacdm_inv_efunc(double z, double Om0, double Ode0, double Ok0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y, double wp,
    double apiv, double wa):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  cdef Odescl = opz**(3. * (1. + wp + apiv * wa)) * exp(-3. * wa * z / opz)
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 + Ode0 * Odescl, -0.5)

######## w0wzCDM
# No relativistic species
def w0wzcdm_inv_efunc_norel(double z, double Om0, double Ode0, double Ok0,
    double w0, double wz):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1. + w0 - wz)) * exp(-3. * wz * z)
  return pow(opz**2 * (opz * Om0 + Ok0) + Ode0 * Odescl, -0.5)

# Massless neutrinos
def w0wzcdm_inv_efunc_nomnu(double z, double Om0, double Ode0, double Ok0,
    double Or0, double w0, double wz):
  cdef double opz = 1.0 + z
  cdef Odescl = opz**(3. * (1. + w0 - wz)) * exp(-3. * wz * z)
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 +
          Ode0 * Odescl, -0.5)

# With massive neutrinos
def w0wzcdm_inv_efunc(double z, double Om0, double Ode0, double Ok0,
    double Ogamma0, double NeffPerNu, int nmasslessnu, list nu_y, double w0,
    double wz):

  cdef double opz = 1.0 + z
  cdef double Or0 = Ogamma0 * (1.0 + nufunc(opz, NeffPerNu, nmasslessnu, nu_y))
  cdef Odescl = opz**(3. * (1. + w0 - wz)) * exp(-3. * wz * z)
  return pow((((opz * Or0 + Om0) * opz) + Ok0) * opz**2 + Ode0 * Odescl, -0.5)

######## Neutrino relative density function
# Scalar equivalent to FLRW.nu_realative_density in core.py
#  Please see that for further discussion.
# This should only be called with massive neutrinos (e.g., nu_y is not empty)
# Briefly, this is just a numerical fitting function to the true relationship,
#  which is too expensive to want to evaluate directly.  The
#  constants which appear are:
#    p = 1.83  -> numerical fitting constant from Komatsu et al.
#  1/p = 0.54644... -> same constant
#    k = 0.3173 -> another fitting constant
#  7/8 (4/11)^(4/3) = 0.2271... -> fermion/boson constant for neutrino
#                                   contribution -- see any cosmology book
#  The Komatsu reference is: Komatsu et al. 2011, ApJS 192, 18
cdef nufunc(double opz, double NeffPerNu, int nmasslessnu, list nu_y):
  cdef int N = len(nu_y)
  cdef double k = 0.3173 / opz
  cdef double rel_mass_sum = nmasslessnu
  cdef unsigned int i
  for i in range(N):
    rel_mass_sum += pow(1.0 + (k * <double>nu_y[i])**1.83, 0.54644808743)
  return 0.22710731766 * NeffPerNu * rel_mass_sum
