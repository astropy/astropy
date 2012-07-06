from __future__ import absolute_import, division, print_function
"""
Core units classes and functions
"""
# requires python 2.6 or later

import re, keyword
import numpy as np
import fractions
import functools


Fraction = fractions.Fraction

_registry = {}

class UnitsException(Exception) :
    pass

class UnitBase(object) :
    """Abstract Base class for units
    
    Should not be used by users directly
    """
    # may use standard library ABC module instead

    def __init__(self) :
        raise ValueError, "Cannot directly initialize abstract base class"
        pass

    def __pow__(self, p) :
        if isinstance(p, tuple) and len(p)==2:
            p = Fraction(p[0],p[1])
        else:
            # allow two possible floating point fractions, all others illegal
            if not int(2*p) == 2*p:
                raise ValueError, "floating values for unit powers must be integers or integers +0.5"               
        return CompositeUnit(1, [self], [p]).simplify()
        

    def __div__(self, m) :
        if isinstance(m, UnitBase) :
            return CompositeUnit(1, [self, m], [1, -1]).simplify()
        else :
            return CompositeUnit(1.0/m, [self], [1]).simplify()

    def __rdiv__(self, m) :
        return CompositeUnit(m, [self], [-1]).simplify()
        
    def __truediv__(self, m) :
        if isinstance(m, UnitBase) :
            return CompositeUnit(1, [self, m], [1, -1]).simplify()
        else :
            return CompositeUnit(1.0/m, [self], [1]).simplify()

    def __rtruediv__(self, m) :
        return CompositeUnit(m, [self], [-1]).simplify()

    def __mul__(self, m) :
        if hasattr(m, "units") :
            return m*self
        elif isinstance(m, UnitBase) :
            return CompositeUnit(1, [self, m], [1,1]).simplify()
        else :
            return CompositeUnit(m, [self], [1]).simplify()

    def __rmul__(self, m) :
        return CompositeUnit(m, [self], [1]).simplify()

    def __repr__(self) :
        return 'Unit("'+str(self)+'")'

    def __eq__(self, other) :
        try:
            return self.convert_to(other,1)==1.
        except UnitsException :
            return False

    def __ne__(self, other) :
        return not (self==other)

    def __lt__(self, other) :
        return self.convert_to(other,1)<1.

    def __gt__(self, other) :
        return self.convert_to(other,1)>1.

    def __le__(self, other) :
        return self.convert_to(other,1)<=1.

    def __ge__(self, other) :
        return self.convert_to(other,1)>=1.

    def __neg__(self) :
        return self*(-1)

    def simplify(self) :
        return self

    def is_dimensionless(self) :
        return False

    def converter_to(self, other):
        """return the conversion function to convert values from to the specified unit
        
        Parameters
        ----------
        other: unit object or string that can be converted to a unit object
        
        Returns
        -------
        A function that normally expects a single argument
        that is a scalar value or an array of values (or anything that may
        be converted to an array). Subclasses may add extra arguments 
        
        Raise
        -----
        UnitException
            If units are inconsistent
        """

        if isinstance(other, str) :
            other = unit(other)
        try :
            scale = (self/other).dimensionless_constant()
            return lambda val: scale * argcondition(val)
        except UnitsException :
            raise UnitsException, "Not convertible"

    def convert_to(self, other, value):
        """return the converted values in the specified unit
        
        Parameters
        ----------
        other: unit object or string that can be converted to a unit object
        value: scalar int or float, or sequence that can be converted to array
            value(s) in the current unit to be converted to the specified unit
        
        Returns
        -------
        Converted value(s). Input value sequences are returned as numpy arrays
        
        Raise
        -----
        UnitException
            If units are inconsistent
        """

        return self.converter_to(other)(value)

    comment  = '''
    def ratio(self, other, **substitutions) :
        """Get the conversion ratio between this Unit and another
        specified unit.

        Keyword arguments, if specified, give numerical substitutions
        for the named unit. This is most useful for specifying values
        for cosmological quantities like 'a' and 'h', but can also
        be used for any IrreducibleUnit.

        >>> unit()"1 Mpc a").ratio("kpc", a=0.25)
        250.0
        >>> unit("1 Mpc").ratio("Msol")
        UnitsException: not convertible
        >>> unit("1 Mpc").ratio("Msol", kg=25.0, m=50.0)
        3.1028701506345152e-08
        """

        if isinstance(other, str) :
            other = unit(other)
        try :
            return (self/other).dimensionless_constant(**substitutions)
        except UnitsException :
            raise UnitsException, "Not convertible"

    def in_units(self, *a, **kw) :
        """Alias for ratio"""
    
        return self.ratio(*a, **kw)
        '''
    
    def irrep(self) :
        """Return a unit object composed of only irreducible units
        
        Parameters
        ----------
        None
        
        Returns
        -------
        New CompositeUnit object containing only irreducible unit objects
        """
        return self

    def _register_unit(self, st) :
        if st in _registry :
            raise UnitsException, "Unit with this name already exists"
        if "**" in st or "^" in st or " " in st :
            # will cause problems for simple string parser in unit() factory
            raise UnitsException, "Unit names cannot contain '**' or '^' or spaces"
        _registry[st]=self

    def __deepcopy__(self, memo) :
        # This may look odd, but the units conversion will be very
        # broken after deep-copying if we don't guarantee that a given
        # physical unit corresponds to only one instance
        return self



class IrreducibleUnit(UnitBase) :
    """Irreducible units all other units of the same kind are defined in terms of
    
    Examples are meters, seconds, kilograms, coulombs, etc. There is only 
    once instance of such a unit per type.
    """
    def __init__(self, st) :
        """"""
        self._st_rep = st
        self._register_unit(st)


    def __str__(self) :
        return self._st_rep

    def latex(self) :
        """Generate latex representation of unit name

        Returns
        -------
        Latex string
        """
        return r"\mathrm{"+self._st_rep+"}"

    def irrep(self) :
        return CompositeUnit(1, [self], [1])


class Unit(UnitBase) :
    def __init__(self, st, represents) :
        self._st_rep = st
        if isinstance(represents, str) :
            represents = unit(represents)
            
        self._represents = represents
        self._register_unit(st)

    def __str__(self) :
        """Return string representation for unit"""
        return self._st_rep

    def latex(self) :
        """Generate latex representation of unit name
        
        Prefactors are converted into exponent notation. Named units by default
        are represented by the string '\mathrm{unit_name}', although this can
        be overriden be overriden by setting unit_name._latex.
        
        Returns
        -------
        Latex string
        """
        if hasattr(self,'_latex') :
            return self._latex
        return r"\mathrm{"+self._st_rep+"}"

    def irrep(self) :
        """Return a unit object composed of only irreducible units
        
        Parameters
        ----------
        None
        
        Returns
        -------
        New CompositeUnit object containing only irreducible unit objects
        """
        return self._represents.irrep()


class CompositeUnit(UnitBase) :
    def __init__(self, scale, bases, powers) :
        """Create a composite unit using expressions of previously defined units.

        Direct use of this function is not recommended. Instead use the
        factory function unit(...).
        """
        
        if scale==1. :
            scale = 1

        self._scale = scale
        self._bases = bases
        self._powers = powers

    def latex(self) :
        """Generate latex representation of unit name
        
        Prefactors are converted into exponent notation. Named units by default
        are represented by the string '\mathrm{unit_name}', although this can
        be overriden by setting unit_name._latex.
        
        Returns
        -------
        Latex string
        """
        
        if self._scale!=1 :
            x = ("%.2e"%self._scale).split('e')
            s = x[0]
            ex = x[1].lstrip('0+')
            if len(ex)>0 and ex[0]=='-':
                ex = '-'+(ex[1:]).lstrip('0')
            if ex!='' : s+=r"\times 10^{"+ex+"}"
        else :
            s = ""

        for b,p in zip(self._bases, self._powers) :
            if s!="" :
                s+=r"\,"+b.latex()
            else :
                s = b.latex()

            if p!=1 :
                s+="^{"
                s+=str(p)
                s+="}"
        return s


    def __str__(self) :
        """Return string representation for unit"""
        s=None
        if len(self._bases)==0 :
            return "%.2e"%self._scale

        if self._scale!=1 :
            s = "%.2e"%self._scale

        for b,p in zip(self._bases, self._powers) :
            if s is not None :
                s+=" "+str(b)
            else :
                s = str(b)

            if p!=1 :
                s+="**"
                if isinstance(p,Fraction) :
                    s+=str(p)
                else :
                    s+=str(p)
        return s

    def _expand(self, expand_to_irrep=False) :
        """Internal routine to expand any pointers to composite units
        into direct pointers to the base units. If expand_to_irrep is
        True, everything is expressed in irreducible units.
        A _gather will normally be necessary to sanitize the unit
        after an _expand."""

        trash = []



        for i,(b,p) in enumerate(zip(self._bases, self._powers)) :
            if isinstance(b,Unit) and expand_to_irrep :
                b = b._represents.irrep()

            if isinstance(b,CompositeUnit) :
                if expand_to_irrep :
                    b = b.irrep()

                trash.append(i)
                self._scale*=b._scale**p
                for b_sub, p_sub in zip(b._bases, b._powers) :
                    self._bases.append(b_sub)
                    self._powers.append(p_sub*p)

        trash.sort()
        for offset,i in enumerate(trash) :
            del self._bases[i-offset]
            del self._powers[i-offset]


    def _gather(self) :
        """Internal routine to gather together powers of the same base
        units, then order the base units by their power (descending)"""

        trash = []
        bases = list(set(self._bases))
        powers = [sum([p for bi,p in zip(self._bases, self._powers)
                       if bi is b]) \
                  for b in bases]

        bp = sorted(filter(lambda x : x[0]!=0,
                           zip(powers, bases)),
                    reverse=True,
                    key=lambda x: x[0])
                    # Py2 only: cmp=lambda x, y: cmp(x[0], y[0]))

        if len(bp)!=0 :
            self._powers, self._bases = map(list,zip(*bp))
        else :
            self._powers, self._bases = [],[]




    def copy(self) :
        """Create a shallow copy 
        
        The returned copy references exactly the same underlying base units, 
        but where the list of those units can be manipulated separately.
        """
        return CompositeUnit(self._scale, self._bases[:], self._powers[:])

    def __copy__(self) :
        """For compatibility with python copy module"""
        return self.copy()

    def simplify(self) :
        self._expand()
        self._gather()
        return self

    def irrep(self) :
        """Return a unit object composed of only irreducible units"""

        x = self.copy()
        x._expand(True)
        x._gather()
        return x

    def is_dimensionless(self) :
        """True if this unit actually translates into a scalar quantity."""
        x = self.irrep()
        if len(x._powers)==0 :
            return True

    def dimensionless_constant(self) :
        """If this unit is dimensionless, return its scalar quantity.

        Direct use of this method is not recommended. It is generally
        better to use the convert_to or converter_to methods instead.
        """
        
        x = self.irrep()
        c = x._scale
        if x._bases:
            raise UnitsException, "Not dimensionless"
        return c

    def _power_of(self, base) :
        if base in self._bases :
            return self._powers[self._bases.index(base)]
        else :
            return 0
            
    bozo = '''
    def dimensional_project(self, basis_units) :
        """Work out how to express the dimensions of this unit relative to the
        specified list of basis units.

        This is used by the framework when making inferences about sensible units to
        use in various situations.
        
        For example, you can represent a length as an energy divided by a force:

           >>> ynit("23 kpc").dimensional_project(["J", "N"])
           array([1, -1], dtype=object)

        However it's not possible to represent a length by energy alone:

           >>> unit("23 kpc").dimensional_project(["J"])
           UnitsException: Basis units do not span dimensions of specified unit

        This function also doesn't know what to do if the result is ambiguous:

           >>> unit("23 kpc").dimensional_project(["J", "N", "kpc"])
           UnitsException: Basis units are not linearly independent
        
        """

        vec_irrep = [Unit(x).irrep() for x in basis_units]
        me_irrep = self.irrep()
        bases = set(me_irrep._bases)
        for vec in vec_irrep :
            bases.update(vec._bases)

        bases = list(bases)

        matrix = np.zeros((len(bases),len(vec_irrep)),dtype=Fraction)

        for base_i, base in enumerate(bases) :
            for vec_i, vec in enumerate(vec_irrep) :
                matrix[base_i,vec_i] = vec._power_of(base)


        # The matrix calculated above describes the transformation M
        # such that v = M.d where d is the sought-after powers of the
        # specified base vectors, and v is the powers in terms of the
        # base units in the list bases.
        #
        # To invert, since M is possibly rectangular, we use the
        # solution to the least-squares problem [minimize (v-M.d)^2]
        # which is d = (M^T M)^(-1) M^T v.
        #
        # If the solution to that does not solve v = M.d, there is no
        # admissable solution to v=M.d, i.e. the supplied base vectors do not span
        # the requires space.
        #
        # If (M^T M) is singular, the vectors are not linearly independent, so any
        # solution would not be unique.



        M_T_M = np.dot(matrix.transpose(),matrix)

        try :
            M_T_M_inv = util.rational_matrix_inv(M_T_M)
        except np.linalg.linalg.LinAlgError :
            raise UnitsException, "Basis units are not linearly independent"

        my_powers = [me_irrep._power_of(base) for base in bases]


        candidate= np.dot(M_T_M_inv, np.dot(matrix.transpose(), my_powers))

        # Because our method involves a loss of information (multiplying
        # by M^T), we could get a spurious solution. Check this is not the case...


        if any(np.dot(matrix, candidate)!=my_powers) :
            # Spurious solution, meaning the base vectors did not span the
            # units required in the first place.
            raise UnitsException, "Basis units do not span dimensions of specified unit"

        return candidate
        '''

def unit(s) :
    """
    Class factory function for units. 
    
    Given a string s, creates a Unit object.
    
    Parameters
    ----------
    s : string
      The string format is:
          [<scale>] [<unit_name>][**<rational_power>] [[<unit_name>] ... ]

    Returns
    -------
    Unit object

    Examples
    --------
    >>> unit("1.e30 kg")
    >>> unit("kpc**2")
    >>> unit("26.2 m s**-1")
    """

    if isinstance(s, UnitBase):
        return s

    x = s.split()
    try:
        scale = float(x[0])
        del x[0]
    except (ValueError, IndexError) :
        scale = 1.0

    units = []
    powers = []


    for com in x :
        if "**" in com or "^" in com :
            s = com.split("**" if "**" in com else "^")
            try :
                u = _registry[s[0]]
            except KeyError :
                raise ValueError, "Unknown unit "+s[0]
            p = Fraction(s[1])
            if p.denominator is 1 :
                p = p.numerator
        else :
            u = _registry[com]
            p = 1

        units.append(u)
        powers.append(p)

    return CompositeUnit(scale, units, powers)

def argcondition(value):
    """validate value is acceptable for conversion purposes
    
    Will convert into an array if not a scalar, and can be converted 
    into an array
    
    Parameters
    ----------
    value: int or float value, or sequence of such values
        that can be converted into an array if not scalar
        
    Returns
    -------
    Scalar value or numpy array
    
    Raises
    ------
    ValueError
        If value is not as expected
    """
    print (type(value))
    if isinstance(value, float) or isinstance(value, int):
        return value
    else:
        try:
            avalue = np.array(value)
            dt = str(avalue.dtype)
            if not (dt.startswith('int') or dt.startswith('float')):
                raise ValueError, "Must be convertable to int or float array"
            return avalue
        except ValueError:
            raise ValueError, \
            "Value not scalar compatible or convertable into a float or integer array"

def takes_arg_in_units(*args, **orig_kwargs) :
    """
    
    Returns a decorator to create a function which auto-converts input
    to given units.

    **Usage:**
    
    .. code-block:: python

        @takes_arg_in_units((2, "Msol"), (1, "kpc"), ("blob", "erg"))
        def my_function(arg0, arg1, arg2, blob=22) :
           print "Arg 2 is",arg2,"Msol"
           print "Arg 1 is",arg1,"kpc"
           print "blob is",blob,"ergs"



    >>> My_function(22, "1.e30 kg", 23, blob="12 J")
    Input 3 is 0.5 Msol
    Input 2 is 23 kpc

    """

    context_arg = orig_kwargs.get('context_arg',None)
    
    kwargs = filter(lambda x: hasattr(x[0],'__len__'), args)
    args = filter(lambda x: not hasattr(x[0], '__len__'), args)
    
    def decorator_fn(x) :
        @functools.wraps
        def wrapper_fn(*fn_args, **fn_kwargs) :
            context = {}
            if context_arg is not None :
                context = fn_args[context_arg].conversion_context()
                
            fn_args = list(fn_args)

            for arg_num, arg_units in args :

                if isinstance(fn_args[arg_num],str) :
                    fn_args[arg_num] = unit(fn_args[arg_num])

                if hasattr(fn_args[arg_num], "in_units") :
                    fn_args[arg_num] = fn_args[arg_num].in_units(arg_units,**context)

            for arg_name, arg_units in kwargs :
                if isinstance(fn_kwargs[arg_name],str) :
                    fn_kwargs[arg_name] = unit(fn_kwargs[arg_name])

                if hasattr(fn_kwargs[arg_name], "in_units") :
                    fn_kwargs[arg_name] = fn_kwargs[arg_name].in_units(arg_units, **context)

            return x(*fn_args, **fn_kwargs)

        return wrapper_fn

    return decorator_fn

for irr_unit_name in ['m','s','kg','K','a','h']:
    globals()[irr_unit_name] = IrreducibleUnit(irr_unit_name)

# Times
minutes = Unit('minutes', 60 * s)
hr = Unit('hr', 3600 * s)

yr = Unit('yr', 3.1556926e7 * s)
kyr = Unit('kyr', 1000 * yr)
Myr = Unit('Myr', 1000 * kyr)
Gyr = Unit('Gyr', 1000 * Myr)

# Frequency

Hz = Unit('Hz', 1 / s)

# Distances
cm = Unit('cm', 0.01 * m)
km = Unit('km', 1000 * m)
au = Unit('au', 1.49598e11 * m)
pc = Unit('pc', 3.08568025e16 * m)
kpc = Unit('kpc', 1000 * pc)
Mpc = Unit('Mpc', 1000 * kpc)
Gpc = Unit('Gpc', 1000 * Mpc)

# Masses
Msol = Unit('Msol', 1.98892e30 * kg)
_registry['Msol']._latex = r'M_{\odot}'
g = Unit('g', 1.0e-3 * kg)
m_p = Unit('m_p', 1.67262158e-27 * kg)
_registry['m_p']._latex = 'm_p'
m_e = Unit('m_e', 9.10938188e-31 * kg)
_registry['m_e']._latex = 'm_e'
# Forces
N = Unit('N', kg * m * s**-2)

# Energies
J = Unit('J', N * m)
erg = Unit('erg', 1.0e-7 * J)
eV = Unit('eV', 1.60217646e-19 * J)
keV = Unit('keV', 1000 * eV)
MeV = Unit('MeV', 1000 * keV)

# Pressures
Pa = Unit('Pa', J * m**-3)
dyn = Unit('dyn', erg * cm**-3)

class EquivalenceUnit(UnitBase):
    """Class to handle equivalence relations between dissimilar units
    
    Must be subclassed
    
    This works by containing (not inheriting) an instance of a specific unit,
    and having a list of equivalent units that are acceptable for conversions.
    The first of these is the one used as the basis for all conversions.
    Should there be more than two equivalent types, the first is the 
    intermediate unit used to convert from the others.
    
    This class is only used to convert between dissimilar units with implicit
    relations. It cannot be used as part of other composite units. Towards
    that end, the functionality of that class for arithmetic operations with
    other unit classes is disabled. (Operations with scalars is permitted.)
    """
    
    def __init__(eunit, name=None):
        raise NotImplementedError, "object initialization must be handled by subclass"
    def __str__(self):
        return str(self.eunit)
    def __repr__(self):
        return "EquivalenceUnit(%s)" % str(self.eunit) # not right
    def _equivalent(self, nunit):
        """which of the internal representations is given unit is compatible with, if any
        
        returns index of the equivalency, negative if none
        """
        for i, eunit in enumerate(self.elist):
            try:
                eunit.converter_to(nunit)
                return i
            except UnitsException:
                pass
        return -1
            
    def convert_to(self, nunit, value, *args, **kwds):
        return self.converter_to(nunit, value, args, kwds)
    
    # doesn't handle extra args or keywords yet!
    def converter_to(self, nunit, *args, **kwds):
        """will need to override if conversion requires associated values"""
        to_index = self._equivalent(self.eunit)
        if isinstance(nunit, self.__class__):
            nunit = nunit.eunit
        from_index = self._equivalent(nunit)
        if from_index < 0:
            raise UnitException, "not permitted to convert to specified type"
        return lambda value: \
            self._from_standard[from_index](nunit, self._to_standard[to_index](value))
     
    def __pow__(self, p) :
        raise NotImplementedError, "powers not permitted for equivalence units"
        
    def __div__(self, m) :
        try:
            factor = float(m)
        except ValueError:
            raise ValueError, "can only divide by scalars"
        return self__class__(self.eunit/float)

    def __rdiv__(self, m) :
        raise NotImplementedError, "inverse not permitted for equivalence units"
        
    def __truediv__(self, m) :
        try:
            factor = float(m)
        except ValueError:
            raise ValueError, "can only divide by scalars"
        return self.__class__(self.eunit/float)

    def __rtruediv__(self, m) :
        raise NotImplementedError, "inverse not permitted for equivalence units"

    def __mul__(self, m) :
        try:
            factor = float(m)
        except ValueError:
            raise ValueError, "can only multiply by scalars"
        return self.__class__(factor*self.eunit)

    def __rmul__(self, m) :
        try:
            factor = float(m)
        except ValueError:
            raise ValueError, "can only multiply by scalars"
        return self.__class__(factor*self.eunit)

# todo: these eventually must come from the Constants module    
c = 2.99792458e8
h = 6.62606957e-34

class SpectralUnit(EquivalenceUnit):
    """Handles spectral wavelength, frequency, and energy equivalences
    
    Allows conversions between wavelength units, frequency units and
    energy units as they relate to light. 
    
    These special units may not be combined with any others. They exist
    purely for conversion between the specified unit representations.
    """
    def __init__(self, eunit):
        self.elist = [m, Hz, J]
        if self._equivalent(eunit) < 0:
            raise ValueError, "unit is not one of the acceptable types"
        self.eunit = eunit
        self._to_standard = [
            lambda value: self.eunit.converter_to(m)(value),
            lambda value: c / self.eunit.converter_to(Hz)(value) ,
            lambda value: h * c / self.eunit.converter_to(J)(value)
        ]
        self._from_standard = [
            lambda nunit, value: m.converter_to(nunit)(value),
            lambda nunit, value: Hz.converter_to(nunit)(c / value),
            lambda nunit, value: J.converter_to(nunit)(h * c / value)
        ]

    def _converter_to_standard(self):
        raise NotImplementedError, "Must be defined by subclass"

    def _converter_from_standard(self, nunit):
           raise NotImplementedError, "Must be defined by subclass"
