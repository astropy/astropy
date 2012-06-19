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
    """Base class for units"""

    def __init__(self) :
        raise ValueError, "Cannot directly initialize abstract base class"
        pass

    def __pow__(self, p) :
        if isinstance(p, tuple) :
            p = Fraction(p[0],p[1])
        return CompositeUnit(1, [self], [p]).simplify()

    def __div__(self, m) :
        if isinstance(m, UnitBase) :
            return CompositeUnit(1, [self, m], [1, -1]).simplify()
        else :
            return CompositeUnit(1.0/m, [self], [1]).simplify()

    def __rdiv__(self, m) :
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
            return self.ratio(other)==1.
        except UnitsException :
            return False

    def __ne__(self, other) :
        return not (self==other)

    def __lt__(self, other) :
        return self.ratio(other)<1.

    def __gt__(self, other) :
        return self.ratio(other)>1.

    def __le__(self, other) :
        return self.ratio(other)<=1.

    def __ge__(self, other) :
        return self.ratio(other)>=1.

    def __neg__(self) :
        return self*(-1)

    def simplify(self) :
        return self

    def is_dimensionless(self) :
        return False

    def converter_to(self, other):
        """return the conversion function between this Unit and another
        specified unit. This function normally expects a single argument
        that is a scalar value or an array of values (or anything that may
        be converted to an array). 
        """

        if isinstance(other, str) :
            other = Unit(other)
        try :
            scale = (self/other).dimensionless_constant()
            return lambda val: scale * argcondition(val)
        except UnitsException :
            raise UnitsException, "Not convertible"

    def convert_to(self, other, value):
        """Convert a value or value array to the specified unit"""
        return self.converter_to(other)(value)

    comment  = '''
    def ratio(self, other, **substitutions) :
        """Get the conversion ratio between this Unit and another
        specified unit.

        Keyword arguments, if specified, give numerical substitutions
        for the named unit. This is most useful for specifying values
        for cosmological quantities like 'a' and 'h', but can also
        be used for any IrreducibleUnit.

        >>> Unit("1 Mpc a").ratio("kpc", a=0.25)
        250.0
        >>> Unit("1 Mpc").ratio("Msol")
        UnitsException: not convertible
        >>> Unit("1 Mpc").ratio("Msol", kg=25.0, m=50.0)
        3.1028701506345152e-08
        """

        if isinstance(other, str) :
            other = Unit(other)
        try :
            return (self/other).dimensionless_constant(**substitutions)
        except UnitsException :
            raise UnitsException, "Not convertible"

    def in_units(self, *a, **kw) :
        """Alias for ratio"""
    
        return self.ratio(*a, **kw)
        '''
    
    def irrep(self) :
        """Return a unit equivalent to this one (may be identical) but
        expressed in terms of the currently defined IrreducibleUnit
        instances."""
        return self

    def _register_unit(self, st) :
        if st in _registry :
            raise UnitsException, "Unit with this name already exists"
        if "**" in st or "^" in st or " " in st :
            # will cause problems for simple string parser in Unit() factory
            raise UnitsException, "Unit names cannot contain '**' or '^' or spaces"
        _registry[st]=self

    def __deepcopy__(self, memo) :
        # This may look odd, but the units conversion will be very
        # broken after deep-copying if we don't guarantee that a given
        # physical unit corresponds to only one instance
        return self



class IrreducibleUnit(UnitBase) :
    def __init__(self, st) :
        self._st_rep = st
        self._register_unit(st)


    def __str__(self) :
        return self._st_rep

    def latex(self) :
        return r"\mathrm{"+self._st_rep+"}"

    def irrep(self) :
        return CompositeUnit(1, [self], [1])


class NamedUnit(UnitBase) :
    def __init__(self, st, represents) :
        self._st_rep = st
        if isinstance(represents, str) :
            represents = Unit(represents)
            
        self._represents = represents
        self._register_unit(st)

    def __str__(self) :
        return self._st_rep

    def latex(self) :
        if hasattr(self,'_latex') :
            return self._latex
        return r"\mathrm{"+self._st_rep+"}"

    def irrep(self) :
        return self._represents.irrep()


class CompositeUnit(UnitBase) :
    def __init__(self, scale, bases, powers) :
        """Initialize a composite unit.

        Direct use of this function is not recommended. Instead use the
        factory function Unit(...)."""
        
        if scale==1. :
            scale = 1

        self._scale = scale
        self._bases = bases
        self._powers = powers

    def latex(self) :
        """Returns a LaTeX representation of this unit.

        Prefactors are converted into exponent notation. Named units by default
        are represented by the string '\mathrm{unit_name}', although this can
        be overriden in the pynbody configuration files or by setting
        unit_name._latex."""
        
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
            if isinstance(b,NamedUnit) and expand_to_irrep :
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
        """Create a copy which is 'shallow' in the sense that it
        references exactly the same underlying base units, but where
        the list of those units can be manipulated separately."""
        return CompositeUnit(self._scale, self._bases[:], self._powers[:])

    def __copy__(self) :
        """For compatibility with python copy module"""
        return self.copy()

    def simplify(self) :
        self._expand()
        self._gather()
        return self

    def irrep(self) :
        """Return a new unit which represents this unit expressed
        solely in terms of IrreducibleUnit bases."""
        x = self.copy()
        x._expand(True)
        x._gather()
        return x

    def is_dimensionless(self) :
        """Returns true if this unit actually translates into a scalar
        quantity."""
        x = self.irrep()
        if len(x._powers)==0 :
            return True

    def dimensionless_constant(self) :
        """If this unit is dimensionless, return its scalar quantity.

        Direct use of this function is not recommended. It is generally
        better to use the ratio function instead.
        
        Provide keyword arguments to set values for named IrreducibleUnits --
        see the ratio function for more information."""
        
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

           >>> Unit("23 kpc").dimensional_project(["J", "N"])
           array([1, -1], dtype=object)

        However it's not possible to represent a length by energy alone:

           >>> Unit("23 kpc").dimensional_project(["J"])
           UnitsException: Basis units do not span dimensions of specified unit

        This function also doesn't know what to do if the result is ambiguous:

           >>> Unit("23 kpc").dimensional_project(["J", "N", "kpc"])
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
        
def Unit(s) :
    """
    Class factory for units. Given a string s, creates
    a Unit object.

    The string format is:
      [<scale>] [<unit_name>][**<rational_power>] [[<unit_name>] ... ]

    for example:

      "1.e30 kg"

      "kpc**2"

      "26.2 m s**-1"
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
                    fn_args[arg_num] = Unit(fn_args[arg_num])

                if hasattr(fn_args[arg_num], "in_units") :
                    fn_args[arg_num] = fn_args[arg_num].in_units(arg_units,**context)

            for arg_name, arg_units in kwargs :
                if isinstance(fn_kwargs[arg_name],str) :
                    fn_kwargs[arg_name] = Unit(fn_kwargs[arg_name])

                if hasattr(fn_kwargs[arg_name], "in_units") :
                    fn_kwargs[arg_name] = fn_kwargs[arg_name].in_units(arg_units, **context)

            return x(*fn_args, **fn_kwargs)

        return wrapper_fn

    return decorator_fn

for irr_unit_name in ['m','s','kg','K','a','h']:
    globals()[irr_unit_name] = IrreducibleUnit(irr_unit_name)

# Times
yr = NamedUnit('yr', 3.1556926e7 * s)
kyr = NamedUnit('kyr', 1000 * yr)
Myr = NamedUnit('Myr', 1000 * kyr)
Gyr = NamedUnit('Gyr', 1000 * Myr)

# Distances
cm = NamedUnit('cm', 0.01 * m)
km = NamedUnit('km', 1000 * m)
au = NamedUnit('au', 1.49598e11 * m)
pc = NamedUnit('pc', 3.08568025e16 * m)
kpc = NamedUnit('kpc', 1000 * pc)
Mpc = NamedUnit('Mpc', 1000 * kpc)
Gpc = NamedUnit('Gpc', 1000 * Mpc)

# Masses
Msol = NamedUnit('Msol', 1.98892e30 * kg)
_registry['Msol']._latex = r'M_{\odot}'
g = NamedUnit('g', 1.0e-3 * kg)
m_p = NamedUnit('m_p', 1.67262158e-27 * kg)
_registry['m_p']._latex = 'm_p'
m_e = NamedUnit('m_e', 9.10938188e-31 * kg)
_registry['m_e']._latex = 'm_e'
# Forces
N = NamedUnit('N', kg * m * s**-2)

# Energies
J = NamedUnit('J', N * m)
erg = NamedUnit('erg', 1.0e-7 * J)
eV = NamedUnit('eV', 1.60217646e-19 * J)
keV = NamedUnit('keV', 1000 * eV)
MeV = NamedUnit('MeV', 1000 * keV)

# Pressures
Pa = NamedUnit('Pa', J * m**-3)
dyn = NamedUnit('dyn', erg * cm**-3)
