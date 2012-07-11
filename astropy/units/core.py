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

class UnitsException(Exception):
    pass

class UnitBase(object):
    """Abstract Base class for units
    
    Should not be used by users directly
    """
    # may use standard library ABC module instead

    def __init__(self):
        raise ValueError, "Cannot directly initialize abstract base class"
        pass

    def __pow__(self, p):
        if isinstance(p, tuple) and len(p)==2:
            p = Fraction(p[0],p[1])
        else:
            # allow two possible floating point fractions, all others illegal
            if not int(2*p) == 2*p:
                raise ValueError, "floating values for unit powers must be integers or integers +0.5"               
        return CompositeUnit(1, [self], [p]).simplify()
        

    def __div__(self, m):
        if isinstance(m, EquivalenceUnit):
            raise TypeError, "cannot combine equivalence unit types with any others"
        if isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1, -1]).simplify()
        else :
            return CompositeUnit(1.0/m, [self], [1]).simplify()

    def __rdiv__(self, m):
        return CompositeUnit(m, [self], [-1]).simplify()
        
    def __truediv__(self, m):
        if isinstance(m, EquivalenceUnit):
            raise TypeError, "cannot combine equivalence unit types with any others"
        if isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1, -1]).simplify()
        else :
            return CompositeUnit(1.0/m, [self], [1]).simplify()

    def __rtruediv__(self, m):
        return CompositeUnit(m, [self], [-1]).simplify()

    def __mul__(self, m):
        if isinstance(m, EquivalenceUnit):
            raise TypeError, "cannot combine equivalence unit types with any others"
        if hasattr(m, "units"):
            return m*self
        elif isinstance(m, UnitBase):
            return CompositeUnit(1, [self, m], [1,1]).simplify()
        else :
            return CompositeUnit(m, [self], [1]).simplify()

    def __rmul__(self, m):
        return CompositeUnit(m, [self], [1]).simplify()

    def __repr__(self):
        return 'unit("'+str(self)+'")'

    def __eq__(self, other):
        try:
            return self.convert_to(other,1)==1.
        except UnitsException :
            return False

    def __ne__(self, other):
        return not (self==other)

    def __lt__(self, other):
        return self.convert_to(other,1)<1.

    def __gt__(self, other):
        return self.convert_to(other,1)>1.

    def __le__(self, other):
        return self.convert_to(other,1)<=1.

    def __ge__(self, other):
        return self.convert_to(other,1)>=1.

    def __neg__(self):
        return self*(-1)

    def simplify(self):
        return self

    def is_dimensionless(self):
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
        UnitsException
            If units are inconsistent
        """

        if isinstance(other, str):
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
    
    def scale_to(self, other):
        """return scale factor for converting current units value to new units value"""
        return self.convert_to(other,1)

    comment  = '''
    def ratio(self, other, **substitutions):
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

        if isinstance(other, str):
            other = unit(other)
        try :
            return (self/other).dimensionless_constant(**substitutions)
        except UnitsException :
            raise UnitsException, "Not convertible"

    def in_units(self, *a, **kw):
        """Alias for ratio"""
    
        return self.ratio(*a, **kw)
        '''
    
    def irrep(self):
        """Return a unit object composed of only irreducible units
        
        Parameters
        ----------
        None
        
        Returns
        -------
        New CompositeUnit object containing only irreducible unit objects
        """
        return self

    def _register_unit(self, var=False):
        if not self._st_rep:
            raise UnitsException, "unit has no string representation"
        namelist = [self._st_rep] + self._aliases
        for st in namelist:
            if st in _registry :
                raise UnitsException, "Unit with this name already exists"
            if "**" in st or "^" in st or " " in st :
                # will cause problems for simple string parser in unit() factory
                raise UnitsException, "Unit names cannot contain '**' or '^' or spaces"
            _registry[st] = self
            if var:
                globals()[st] = self

    def __deepcopy__(self, memo):
        # This may look odd, but the units conversion will be very
        # broken after deep-copying if we don't guarantee that a given
        # physical unit corresponds to only one instance
        return self



class IrreducibleUnit(UnitBase):
    """Irreducible units all other units of the same kind are defined in terms of
    
    Examples are meters, seconds, kilograms, coulombs, etc. There is only 
    once instance of such a unit per type.
    """
    def __init__(self, st):
        """"""
        if isinstance(st, str):
            self._aliases = []
            self._st_rep = st
        else:
            if len(st) == 0:
                raise ValueError, "alias list must have at least one entry"
            try:
                self._st_rep = st[0]
                self._aliases = st[1:]
            except IndexError:
                raise ValueError, "name argument must be a string or a list of strings"
        self._register_unit()


    def __str__(self):
        return self._st_rep

    def latex(self):
        """Generate latex representation of unit name

        Returns
        -------
        Latex string
        """
        return r"\mathrm{"+self._st_rep+"}"

    def irrep(self):
        return CompositeUnit(1, [self], [1])


class Unit(UnitBase):
    def __init__(self, st, represents, var=False):
        """Create a named unit. 
        
        if var is True, create variables in namespace for each alias
        """
        if isinstance(st, str):
            self._aliases = []
            self._st_rep = st
        else:
            if len(st) == 0:
                raise ValueError, "alias list must have at least one entry"
            try:
                self._st_rep = st[0]
                self._aliases = st[1:]
            except IndexError:
                raise ValueError, "name argument must be a string or a list of strings"
        if isinstance(represents, str):
            represents = unit(represents)           
        self._represents = represents
        self._register_unit(var)

    def __str__(self):
        """Return string representation for unit"""
        return self._st_rep

    def latex(self):
        """Generate latex representation of unit name
        
        Prefactors are converted into exponent notation. Named units by default
        are represented by the string '\mathrm{unit_name}', although this can
        be overriden be overriden by setting unit_name._latex.
        
        Returns
        -------
        Latex string
        """
        if hasattr(self,'_latex'):
            return self._latex
        return r"\mathrm{"+self._st_rep+"}"

    def irrep(self):
        """Return a unit object composed of only irreducible units
        
        Parameters
        ----------
        None
        
        Returns
        -------
        New CompositeUnit object containing only irreducible unit objects
        """
        return self._represents.irrep()


class CompositeUnit(UnitBase):
    def __init__(self, scale, bases, powers):
        """Create a composite unit using expressions of previously defined units.

        Direct use of this function is not recommended. Instead use the
        factory function unit(...).
        """
        
        if scale==1. :
            scale = 1

        self._scale = scale
        self._bases = bases
        self._powers = powers

    def latex(self):
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

        for b,p in zip(self._bases, self._powers):
            if s!="" :
                s+=r"\,"+b.latex()
            else :
                s = b.latex()

            if p!=1 :
                s+="^{"
                s+=str(p)
                s+="}"
        return s


    def __str__(self):
        """Return string representation for unit"""
        s=None
        if len(self._bases)==0 :
            return "%.2e"%self._scale

        if self._scale!=1 :
            s = "%.2e"%self._scale

        for b,p in zip(self._bases, self._powers):
            if s is not None :
                s+=" "+str(b)
            else :
                s = str(b)

            if p!=1 :
                s+="**"
                if isinstance(p,Fraction):
                    s+=str(p)
                else :
                    s+=str(p)
        return s

    def _expand(self, expand_to_irrep=False):
        """Internal routine to expand any pointers to composite units
        into direct pointers to the base units. If expand_to_irrep is
        True, everything is expressed in irreducible units.
        A _gather will normally be necessary to sanitize the unit
        after an _expand."""

        trash = []



        for i,(b,p) in enumerate(zip(self._bases, self._powers)):
            if isinstance(b,Unit) and expand_to_irrep :
                b = b._represents.irrep()

            if isinstance(b,CompositeUnit):
                if expand_to_irrep :
                    b = b.irrep()

                trash.append(i)
                self._scale*=b._scale**p
                for b_sub, p_sub in zip(b._bases, b._powers):
                    self._bases.append(b_sub)
                    self._powers.append(p_sub*p)

        trash.sort()
        for offset,i in enumerate(trash):
            del self._bases[i-offset]
            del self._powers[i-offset]


    def _gather(self):
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




    def copy(self):
        """Create a shallow copy 
        
        The returned copy references exactly the same underlying base units, 
        but where the list of those units can be manipulated separately.
        """
        return CompositeUnit(self._scale, self._bases[:], self._powers[:])

    def __copy__(self):
        """For compatibility with python copy module"""
        return self.copy()

    def simplify(self):
        self._expand()
        self._gather()
        return self

    def irrep(self):
        """Return a unit object composed of only irreducible units"""

        x = self.copy()
        x._expand(True)
        x._gather()
        return x

    def is_dimensionless(self):
        """True if this unit actually translates into a scalar quantity."""
        x = self.irrep()
        if len(x._powers)==0 :
            return True

    def dimensionless_constant(self):
        """If this unit is dimensionless, return its scalar quantity.

        Direct use of this method is not recommended. It is generally
        better to use the convert_to or converter_to methods instead.
        """
        
        x = self.irrep()
        c = x._scale
        if x._bases:
            raise UnitsException, "Not dimensionless"
        return c

    def _power_of(self, base):
        if base in self._bases :
            return self._powers[self._bases.index(base)]
        else :
            return 0
            
    bozo = '''
    def dimensional_project(self, basis_units):
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

        for base_i, base in enumerate(bases):
            for vec_i, vec in enumerate(vec_irrep):
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


        if any(np.dot(matrix, candidate)!=my_powers):
            # Spurious solution, meaning the base vectors did not span the
            # units required in the first place.
            raise UnitsException, "Basis units do not span dimensions of specified unit"

        return candidate
        '''

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

    def __init__(self, st, represents):
        raise NotImplementedError, "object initialization must be handled by subclass"
    def __str__(self):
        return self._st_rep
    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, str(self.eunit))
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
        return self.converter_to(nunit, value, *args, **kwds)(value, *args, **kwds)

    def converter_to(self, nunit, *args, **kwds):
        """will need to override if conversion requires associated values"""
        to_index = self._equivalent(self.eunit)
        if isinstance(nunit, self.__class__):
            nunit = nunit.eunit
        from_index = self._equivalent(nunit)
        if from_index < 0:
            raise UnitsException, "not permitted to convert to specified type"
        return lambda value, *args, **kwds: \
            self._from_standard[from_index](nunit, self._to_standard[to_index]
                          (value, *args, **kwds), *args, **kwds)

    def __pow__(self, p):
        raise NotImplementedError, "powers not permitted for equivalence units"

    def __div__(self, m):
        try:
            factor = float(m)
        except TypeError:
            raise TypeError, "can only divide by scalars"
        return self__class__(self.eunit/float)

    def __rdiv__(self, m):
        raise NotImplementedError, "inverse not permitted for equivalence units"

    def __truediv__(self, m):
        try:
            factor = float(m)
        except TypeError:
            raise TypeError, "can only divide by scalars"
        return self.__class__(self.eunit/float)

    def __rtruediv__(self, m):
        raise NotImplementedError, "inverse not permitted for equivalence units"

    def __mul__(self, m):
        try:
            factor = float(m)
        except TypeError:
            raise TypeError, "can only multiply by scalars"
        return self.__class__(factor*self.eunit)

    def __rmul__(self, m):
        try:
            factor = float(m)
        except TypeError:
            raise TypeError, "can only multiply by scalars"
        return self.__class__(factor*self.eunit)


def unit(s):
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
    except (ValueError, IndexError):
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
        if len(units) > 1:
            for u in units:
                if isinstance(u, EquivalenceUnit):
                    raise TypeError, "cannot combine equivalence units with any others"
    if len(units) == 1 and isinstance(units[0],EquivalenceUnit):
        return units[0]
    else:
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

def takes_arg_in_units(*args, **orig_kwargs):
    """
    
    Returns a decorator to create a function which auto-converts input
    to given units.

    **Usage:**
    
    .. code-block:: python

        @takes_arg_in_units((2, "Msol"), (1, "kpc"), ("blob", "erg"))
        def my_function(arg0, arg1, arg2, blob=22):
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
    
    def decorator_fn(x):
        @functools.wraps
        def wrapper_fn(*fn_args, **fn_kwargs):
            context = {}
            if context_arg is not None :
                context = fn_args[context_arg].conversion_context()
                
            fn_args = list(fn_args)

            for arg_num, arg_units in args :

                if isinstance(fn_args[arg_num],str):
                    fn_args[arg_num] = unit(fn_args[arg_num])

                if hasattr(fn_args[arg_num], "in_units"):
                    fn_args[arg_num] = fn_args[arg_num].in_units(arg_units,**context)

            for arg_name, arg_units in kwargs :
                if isinstance(fn_kwargs[arg_name],str):
                    fn_kwargs[arg_name] = unit(fn_kwargs[arg_name])

                if hasattr(fn_kwargs[arg_name], "in_units"):
                    fn_kwargs[arg_name] = fn_kwargs[arg_name].in_units(arg_units, **context)

            return x(*fn_args, **fn_kwargs)

        return wrapper_fn

    return decorator_fn

for irr_unit_name in [['m','meter'],['s','second'],['kg','kilogram'],['K','Kelvin'],'a','h','coulomb']:
    if isinstance(irr_unit_name, str):
        irrunit = IrreducibleUnit(irr_unit_name)
        globals()[irr_unit_name] = irrunit
        _registry[irrname] = irrunit
    else:
        irrunit = IrreducibleUnit(irr_unit_name[0])
        for irrname in irr_unit_name:
            globals()[irrname] = irrunit
            _registry[irrname] = irrunit

def list_like(u):
    """List all the units that are the same type as the specified unit.
    
    Any aliases are noted as such.
    Equivalance units are excluded currently"""
    likedict = {}
    irraliases = []
    for ukey in _registry:
        try:
            tunit = _registry[ukey]
            tunit.converter_to(u)
            if tunit._st_rep == ukey and not isinstance(tunit, EquivalenceUnit):
                likedict[ukey] = tunit
            if isinstance(tunit, IrreducibleUnit) and tunit._st_rep != ukey:
                irraliases.append(ukey)
        except UnitsException:
            pass
    if len(likedict) == 0:
        print("No similar units found")
    else:
        print("Primary name | Unit definition | Aliases" )
        for ukey in likedict:
            uv = likedict[ukey]
            if isinstance(uv, IrreducibleUnit):
                udef = "irreducible"
                ualiases =irraliases
            else:
                udef = uv._represents
                ualiases = uv._aliases
            print("%-15s %-15s %s" % (uv._st_rep,udef,ualiases)) 
                
               

# Times,var=True)
Unit('minutes', 60 * s,var=True)
Unit(['ms','millisecond'], 0.001 * s,var=True)
Unit(['us','microsecond'], 0.001 * ms,var=True)
Unit(['ps','picosecond'], 0.001 * us,var=True)
Unit(['fs','femptosecond'], 0.001 * ps,var=True)

Unit(['hr','hour'], 3600 * s,var=True)
Unit(['day'],24*hr,var=True)
Unit(['sday','sideral_day'],86164.09053*s,var=True)
Unit(['wk','week'], 7*day,var=True)
Unit(['fortnight'], 2*wk,var=True)
Unit(['yr','year'], 3.1556926e7 * s,var=True)
Unit(['Kyr','Kyear','kiloyear'], 1000 * yr,var=True)
Unit(['Myr','Myear','megayear'], 1000 * Kyr,var=True)
Unit(['Gyr','Gyear','gigayear'], 1000 * Myr,var=True)

# Frequency

Unit(['Hz','Hertz','hertz'], 1 / s,var=True)
Unit(['KHz','KHertz','kilohertz'], 1000 * Hz,var=True)
Unit(['MHz','MHertz','megahertz'], 1000 * KHz,var=True)
Unit(['GHz','GHertz','gigahertz'], 1000 * MHz,var=True)

# Distances

Unit(['mm','millimeter'], 0.001 * m,var=True)
Unit(['cm','centimeter'], 0.01 * m,var=True)
Unit(['km','kilometer'], 1000 * m,var=True)
Unit(['au','astronomical_unit'], 1.49598e11 * m,var=True)
Unit(['pc','parsec'], 3.08568025e16 * m,var=True)
Unit(['Kpc','Kparsec','kiloparsec'], 1000 * pc,var=True)
Unit(['Mpc','Mparsec','megaparsec'], 1000 * Kpc,var=True)
Unit(['Gpc','Gparsec','gigaparsec'], 1000 * Mpc,var=True)
Unit(['um','micron'], 0.001 * mm,var=True)
Unit(['nm','nanometer'], 0.001 * um,var=True)
Unit(['A','Angstrom','angstrom'], 0.1 * nm,var=True)
Unit(['pm','picometer'], 0.001 * nm,var=True)
Unit(['inch'], 2.54 * cm,var=True)
Unit(['ft','foot'], 12*inch,var=True)
Unit(['yd','yard'], 3*ft,var=True)
Unit(['mi','mile'], 5280*ft,var=True)

# Masses
Unit('Msol', 1.98892e30 * kg,var=True)
_registry['Msol']._latex = r'M_{\odot}'
Unit(['g','gram'], 1.0e-3 * kg,var=True)
Unit('m_p', 1.67262158e-27 * kg,var=True)
_registry['m_p']._latex = 'm_p'
Unit('m_e', 9.10938188e-31 * kg,var=True)
_registry['m_e']._latex = 'm_e'
Unit(['oz','ounce'],28.349523125*g,var=True) # well, force actually, but who uses it that way?
Unit(['lb','pound'],16*oz,var=True)
Unit(['ton'],2000*lb,var=True)

# Charge
Unit(['e'],coulomb/6.24150965e18,var=True)
# Current
Unit(['amp','ampere'],coulomb/s,var=True)
Unit(['ma','milliamp','milliampere'],0.001*amp,var=True)
# Forces
Unit(['N','Newton','newton'], kg * m * s**-2,var=True)
# Areas
Unit(['barn'],10**-28*m**2,var=True)
Unit(['acre'],43560*ft**2,var=True)

# Volumes
Unit(['l','liter'],1000*cm**3,var=True)
Unit(['ml','milliliter','cc'],cm**3,var=True)
Unit(['gallon'],liter/0.264172052,var=True)
Unit(['quart'],gallon/4,var=True)
Unit(['pint'],quart/2,var=True)
Unit(['cup'],pint/2,var=True)
Unit(['foz','fluid_oz','fluid_ounce'],cup/8,var=True)
Unit(['tbsp','tablespoon'],foz/2,var=True)
Unit(['tsp','teaspoon'],tbsp/3,var=True)
# Energies
Unit(['J','Joule','joule'], N * m,var=True)
Unit('erg', 1.0e-7 * J,var=True)
Unit(['eV','electron_volt'], 1.60217646e-19 * J,var=True)
Unit(['KeV','Kelectron_volt','kilo_electron_volt'], 1000 * eV,var=True)
Unit(['MeV','Melectron_volt','mega_electron_volt'], 1000 * KeV,var=True)
Unit(['BTU','btu'],1.05505585e10*erg,var=True)
Unit(['cal','calorie'],41840000*erg,var=True)
Unit(['kcal','Cal','Calorie','kilocal','kilocalorie'],1000*cal,var=True)
# Power
Unit(['W','Watt','watt'],J/s,var=True)
Unit(['hp','horsepower'],W/0.00134102209,var=True)
# Pressures
Unit('Pa', J * m**-3,var=True)
Unit('dyn', erg * cm**-3,var=True)

# Spectral density

Unit(['flam','flambda'],erg/angstrom/cm**2/s,var=True)
Unit(['fnu'],erg/Hz/cm**2/s,var=True)
Unit(['Jy','Jansky','jansky'],10**-23*fnu,var=True)
Unit(['mJy','mJansky','millijansky'],0.001*Jy,var=True)
Unit(['suJy','uJansky','microjansky'],0.001*mJy,var=True)

# todo: these eventually must come from the Constants module    
c = 2.99792458e8
c_Aps = c * 10**10
h = 6.62606957e-34

class SpectralUnit(EquivalenceUnit):
    """Handles spectral wavelength, frequency, and energy equivalences
    
    Allows conversions between wavelength units, frequency units and
    energy units as they relate to light. 
    
    These special units may not be combined with any others. They exist
    purely for conversion between the specified unit representations.
    
    Note to python novices: use of "lambda" in the code has nothing
    to do with wavelength; see the python lambda statement
    """
    def __init__(self, st, represents, var=False):
        if isinstance(st, str):
            self._aliases = []
            self._st_rep = st
        else:
            if len(st) == 0:
                raise ValueError, "alias list must have at least one entry"
            try:
                self._st_rep = st[0]
                self._aliases = st[1:]
            except IndexError:
                raise ValueError, "name argument must be a string or a list of strings"
        if isinstance(represents, str):
            represents = unit(represents)
        self._represents = represents
        self._register_unit(var)
        self.elist = [m, Hz, J]
        if self._equivalent(represents) < 0:
            raise ValueError, "unit is not one of the acceptable types"
        self.eunit = represents
        self._to_standard = [
            lambda value: self.eunit.converter_to(m)(value),
            lambda value: c / self.eunit.converter_to(Hz)(value),
            lambda value: h * c / self.eunit.converter_to(J)(value)
        ]
        self._from_standard = [
            lambda nunit, value: m.converter_to(nunit)(value),
            lambda nunit, value: Hz.converter_to(nunit)(c / value),
            lambda nunit, value: J.converter_to(nunit)(h * c / value)
        ]

SpectralUnit(['sp_A','sp_Angstrom','sp_angstrom'],A,var=True)
SpectralUnit(['sp_nm','sp_nanometer'],nm,var=True)
SpectralUnit(['sp_um','sp_micron'],um,var=True)
SpectralUnit(['sp_m','sp_meter'],m,var=True)
SpectralUnit(['sp_mm','sp_millimeter'],mm,var=True)
SpectralUnit(['sp_cm','sp_centimeter'],cm,var=True)
SpectralUnit(['sp_Hz','sp_Hertz','sp_hertz'],Hz,var=True)
SpectralUnit(['sp_KHz','sp_KHertz','sp_kilohertz'],KHz,var=True)
SpectralUnit(['sp_MHz','sp_MHertz','sp_megahertz'],MHz,var=True)
SpectralUnit(['sp_GHz','sp_GHertz','sp_gigahertz'],GHz,var=True)
SpectralUnit(['sp_eV','sp_electron_volt'],eV,var=True)
SpectralUnit(['sp_KeV','sp_Kelectron_volt','sp_kilo_electron_volt'],KeV,var=True)
SpectralUnit(['sp_MeV','sp_Melectron_volt','sp_mega_electron_volt'],MeV,var=True)

class SpectralDensityUnit(EquivalenceUnit):
    """Handles spectral density regards wavelength and frequency

    These special units may not be combined with any others. They exist
    purely for conversion between the specified unit representations.
    
    Note to python novices: use of "lambda" in the code has nothing
    to do with wavelength; see the python lambda statement
    """
    def __init__(self, st, represents, var=False):
        if isinstance(st, str):
            self._aliases = []
            self._st_rep = st
        else:
            if len(st) == 0:
                raise ValueError, "alias list must have at least one entry"
            try:
                self._st_rep = st[0]
                self._aliases = st[1:]
            except IndexError:
                raise ValueError, "name argument must be a string or a list of strings"
        if isinstance(represents, str):
            represents = unit(represents)
        self._represents = represents
        self._register_unit(var)
        self.elist = [erg/angstrom/cm**2/s, erg/Hz/cm**2/s]
        if self._equivalent(represents) < 0:
            raise ValueError, "unit is not one of the acceptable types"
        self.eunit = represents
        self._to_standard = [
            lambda value, sunit, sp: (
                self.eunit.converter_to(self.elist[0])(value)
                ),
            lambda value, sunit, sp: (
                c_Aps * self.eunit.converter_to(self.elist[1])(value)
                  / (sunit.converter_to(sp_A)(sp))**2
                )
        ]
        self._from_standard = [
            lambda nunit, value, sunit, sp: (
                self.elist[0].converter_to(nunit)(value)
                ),
            lambda nunit, value, sunit, sp: (
                self.elist[1].converter_to(nunit)(value) * (sunit.converter_to(sp_A)(sp)**2 / c_Aps)
                ),
        ]
        
SpectralDensityUnit(['sd_A','sd_Angstrom','sd_angstrom','sd_flam','sd_flambda'],flambda,var=True)
SpectralDensityUnit(['sd_nm','sd_nanometer'],0.1*flambda,var=True)
SpectralDensityUnit(['sd_um','sd_micron'],0.0001*flambda,var=True)
SpectralDensityUnit(['sd_Hz','sd_Hertz','sd_hertz','sd_fnu'],fnu,var=True)
SpectralDensityUnit(['sd_Jy','sd_Jansky','sd_jansky'],Jy,var=True)
SpectralDensityUnit(['sd_mJy','sd_mJansky','sd_milli_jansky'],0.001*Jy,var=True)
SpectralDensityUnit(['sd_uJy','sd_uJansky','sd_micro_jansky'],0.001*mJy,var=True)

