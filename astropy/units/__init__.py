"""
Module for handling physical units
"""

from __future__ import absolute_import, division, print_function

c_constant = 299792458. # m/s # xxx eventually use from constants
h_constant = 6.62606957e-34 # xxx eventually use from constants

class Unit(object):
    """
    Abstract base class for fundamental units.
    
    Fundamental units are single case of Length, Time, Mass, Charge, etc.
    
    All instances should support the following attributes:
    
    .name
    .scale
    .intrinsic_unit # corresponds to scale = 1.0
    """
    def __init__(self, name='', scale=1.0):
        """Standard initialization"""
        self.name = name
        self.scale = scale
        self.compound_convert = True

    def copy(self):
        return self.__class__(name=self.name,scale=self.scale)

    def isconsistent(self, newunit):
        """Are units compatible?"""
        if not isinstance(newunit, CompoundUnit) and \
            self.__class__ != newunit.__class__:
            return False
        else:
            return True
                
    def convert(self, newunit):
        """Returns a callable object for applying conversions
        
        Only works for fundamental units. Is overridden by CompoundUnit class
        """
        # first check to see if units are compatible
        if self.isconsistent(newunit):
            return lambda x: x * self.scale/newunit.scale
        else:
            raise ValueError("Incompatible unit")
            
    def _factor(self, value):
        try:
            factor = float(value)
        except ValueError:
            raise ValueError("can only multiply by float-compatible numbers")
        return factor
    
    def _combine(self, value, op=None):
        if isinstance(value, Unit):
            unit = CompoundUnit(self, value, op=op)
            if len(unit.num_units) == 1 and len(unit.denom_units) == 0:
                # result is a basic unit so just return a scaled version 
                tunit = unit.num_units[0].copy()
                tunit.scale *= unit.scale
                return tunit
            else:
                return unit
        else:
            if op=="multiply":
                newinst = self.copy()
                newinst.scale *= self._factor(value)
                return newinst
            elif op=="divide":
                newinst = self.copy()
                newinst.scale /= self._factor(value)
                return newinst
            elif op=="rdivide":
                tunit = self.copy()
                tunit.scale /= self._factor(value)
                return CompoundUnit(None, tunit, op='divide')
            elif op=="power":
                pass
            else:
                raise ValueError("unrecognized operation on units")
                   
    def __mul__(self, rhs):
        return self._combine(rhs, op="multiply")
        
    def __rmul__(self, lhs):
        return self._combine(lhs, op="multiply")
        
    def __div__(self, rhs):
        return self._combine(rhs, op="divide")

    def __truediv__(self, rhs):
        return self._combine(rhs, op="divide")

    def __rdiv__(self, lhs):
        return self._combine(lhs, op="rdivide")

    def __rtruediv__(self, lhs):
        return self._combine(lhs, op="rdivide")

        
    def __pow__(self, rhs):
        if type(rhs) != type(0):
            raise ValueError("power must be an integer")
        if rhs == -1:
            return 1./self
        elif rhs < 0: # use * and / recursively for now
            return (1./self)*self**(rhs+1)
        elif rhs == 2:
            return self*self
        elif rhs >= 2:
            return self*self**(rhs-1)
        else:
            return self
            
    def __eq__(self, other):
        return ((self.__class__ == other.__class__) and
                (self.name == other.name) and
                (self.scale == other.scale) and
                (self.compound_convert == other.compound_convert))
        
    def __add__(self, name):
        if type(name) != type(''):
            raise ValueError("items added to units must be strings")
        self.name = name
        return self
        
    def __str__(self):
        if not self.name:
            name = "%f * %s" % (self.scale, self.intrinsic)
        else:
            name = self.name
        return "Units: " + name
        
    def __repr__(self):
        return "%s(scale=%e, name='%s')" % \
            (self.__class__.__name__,self.scale,str(self.name))
        
class CompoundUnit(Unit):
    """Used when the unit involves more than one base unit class
    
    Only recognized operations are multiply, divide, power
    """
    def __init__(self, unit1, unit2=None, op=None, scale=1.0, name=""):
        Unit.__init__(self, name=name)
        if isinstance(unit1, CompoundUnit):
            self.num_units = unit1.num_units[:]
            self.denom_units = unit1.denom_units[:]
            self.scale = unit1.scale
        else:
            if isinstance(unit1,Unit):
                self.num_units = [unit1]
            else:
                self.num_units = []
            self.denom_units = []
        self.scale *= scale
        if unit2 is None:
            unit2_num = []
            unit2_denom =[]
            unit2_scale = 1.
        elif not isinstance(unit2, CompoundUnit):
            unit2_num = [unit2]
            unit2_denom = []
            unit2_scale = 1.
        else:
            unit2_num = unit2.num_units[:]
            unit2_denom = unit2.denom_units[:]
            unit2_scale = unit2.scale
        if op == "multiply":
            self.num_units += unit2_num
            self.denom_units += unit2_denom
            self.scale *= unit2_scale
        elif op == "divide":
            self.num_units += unit2_denom
            self.denom_units += unit2_num
            self.scale /= unit2_scale

        # now order both numerator and denominator in standard ordering
        self.num_units.sort(key=lambda x: x.__class__.__name__)
        self.denom_units.sort(key=lambda x: x.__class__.__name__)
        self._simplify() # "cancel" out identical types of units in numerator
                         # and denominator

    def __eq__(self, other):
        return  ((self.__class__ == other.__class__) and
                  (self.name == other.name) and
                  (self.scale == other.scale) and
                  (self.compound_convert == other.compound_convert) and
                  (self.num_units == other.num_units) and
                  (self.denom_units == other.denom_units))

    def copy(self):
        return self.__class__(self, scale=1., op='multiply', name=self.name)

    def _simplify(self):
        """remove basic types that appear in both numerator and denominator"""
        
        old_num = self.num_units[:]
        old_denom = self.denom_units[:]
        new_num = []
        new_denom = []
        # iterate through both lists by comparing types
        scale = 1.0
        while old_num and old_denom:
            if old_num[0].__class__ == old_denom[0].__class__:
                scale *= 1./old_denom[0].convert(old_num[0])(1)
                del old_num[0]
                del old_denom[0]
            elif (old_num[0].__class__.__name__ < 
                    old_denom[0].__class__.__name__):
                new_num.append(old_num[0])
                del old_num[0]
            else:
                new_denom.append(old_denom[0])
                del old_denom[0]
        new_num += old_num
        new_denom += old_denom
        self.num_units = new_num
        self.denom_units = new_denom
        self.scale *= scale
        
    def _are_matched_lists(self, oldlist, newlist):
        """check two lists of fundamental units are of matching types"""
        for item1, item2 in zip(oldlist, newlist):
            if item1.__class__ != item2.__class__:
                return False
        return True

    def _netscale(self, oldlist, newlist):
        scale = 1.0
        for item1, item2 in zip(oldlist, newlist):
            if item1.compound_convert:
                scale *= item1.convert(item2)(1.)
            else:
                raise ValueError(
                     "Contains a component unit not permitted in compound conversions: %s"
                      % item1.__class__.__name__)
        return scale

    def isconsistent(self, newunit):
        """is the new unit consistent with this object?"""
        if not isinstance(newunit, CompoundUnit) or \
            len(self.num_units) != len(newunit.num_units) or \
            len(self.denom_units) != len(newunit.denom_units):
            return False
        if not self._are_matched_lists(self.num_units,newunit.num_units) or \
            not self._are_matched_lists(self.denom_units,newunit.denom_units):
            return False
        return True
    
    def convert(self, newunit):
        """generalized to compound units"""
        # check for consistency
        if not self.isconsistent(newunit):
            raise ValueError("new unit is inconsistent")
        # determine net scale factor
        factor = (self._netscale(self.num_units,newunit.num_units) /
                  self._netscale(self.denom_units,newunit.denom_units))
        return lambda x: x * factor

    def _build_string(self, sfunc):
        """used for __str__ and __repr__"""
        lnum = len(self.num_units)
        ldenom = len(self.denom_units)
        if self.scale == 1.0:
            ustr = ""
        else:
            ustr = "%e " % self.scale
            if lnum:
                ustr += " * "
        if lnum:
            ustr += " %s " % sfunc(self.num_units[0]) + \
                (lnum-1) * "* %s" % \
                tuple([sfunc(item) for item in self.num_units[1:]])
        else:
            if not ustr:
                ustr = '1 '
        if ldenom:
            udstr = " %s " % sfunc(self.denom_units[0]) + \
                (ldenom-1) * "* %s" % \
                tuple([sfunc(item) for item in self.denom_units[1:]])
        if ldenom == 1:
            ustr += '/ ' + udstr
        elif ldenom > 1:
            ustr += '/ (' + udstr + ')'
        return ustr
       
    def __str__(self):
        return self._build_string(str)
        
    def __repr__(self):
        return self._build_string(repr)

class Length(Unit):
    """
    Standard class for all lengths
    
    Intrinsic length is meters
    """
    def __init__(self, name="", scale=1.0):
        Unit.__init__(self, name, scale)
        self.intrinsic = "meter"
        
meter = m = Length("meter", 1.0)
cm = centimeter = 0.01 * meter + "centimeter"
mm = millimeter = 0.001 * meter + "millimeter"
micron = 10**-6 * meter + "micron"
nm = nanometer = 10**-9 * meter + "nanometer"
angstrom = 10**-10 * meter + "angstrom"

class Mass(Unit):
    """
    Standard class for all masses
    
    Intrinsic mass is kg
    """
    def __init__(self, name="", scale= 1.0):
        Unit.__init__(self, name, scale)
        self.intrinsic = "kilogram"

kilogram = kg = Mass("kilogram",1.0)
gram = 0.001 * kilogram + "gram"

class Time(Unit):
    """
    Standard class for all times
    
    Intrinsic time is seconds
    """
    def __init__(self, name="", scale= 1.0):
        Unit.__init__(self, name, scale)
        self.intrinsic = "second"
        
second = s = Time("second", 1.0)
millisecond = 0.001 * second + "millisecond"
microsecond = 10**-6 * second + "microsecond"
nanosecond = 10**-9 * second + "nanosecond"
picosecond = 10**-12 * second + "picosecond"
femtosecond = 10**-15 * second + "femptosecond"
minute = 60 * second + "minute"
hour = 60 * minute + "hour"
day = 24 * hour + "day"
week = 7 * day + "week"


hertz = hz = 1./second + "hertz"
joule = kg * meter**2 / second**2

spec_map = {'hertz':'frequency',
            'micron':'wavelength',
            'ev':'energy'}

class Spectral(Unit):
    """
    Standard class for all measurements for wavelength/frequency
    
    Handles implicit conversion between length, hz, and energy
    representations.
    """
    def __init__(self, name="", scale= 1.0, flavor=None):
        """This unit comes in 3 flavors wavelength, frequency, and energy
        
        Each of these has a different intrinsic unit depending 
        on how it is initialized
        """
        Unit.__init__(self, name, scale)
        self.compound_convert = False
        if flavor == 'frequency':
            self.intrinsic = 'hertz'
        elif flavor == 'wavelength':
            self.intrinsic = 'micron'
        elif flavor == 'energy':
            self.intrinsic = 'ev'
        else:
            raise ValueError("spectral flavor must be specified")
    
    def flavor(self):
        """Return type of unit"""
        return spec_map[self.instrinsic]  
                      
    def convert(self, newunit):
        """Deal with unit equivalences"""
        if isinstance(newunit, Spectral):
            flavor == self.flavor()
        else:
            if newunit.isconsistent(joule):
                flavor == 'energy'
            elif newunit.isconsisten(meter):
                flavor == 'wavelength'
            elif newunit.isconsistent(hz):
                flavor == 'frequency'
            else:
                raise ValueError(
                  'can only convert to units consistent with'+
                  ' length, frequency or energy')
        