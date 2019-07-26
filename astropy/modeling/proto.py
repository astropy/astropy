'''
Prototyping alternative to use of descriptors for parameters

Need to implement alternative of some descriptor properties.
In particular:

1) be able to define parameters at the class level that have a
name that gets associated with the name of the class attribute.

2) That attribute acts as a  property in that access to the
referenced object is restricted.
'''

class Parm:

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        try:
            self._value = float(value)
        except ValueError:
            raise ValueError("value must be float-compatible")

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    def __init__(self, value=None, default=None, name=None):
        if value is None and default is None:
            raise ValueError("must supply value or default")
        if value is not None:
            self.value = value
        else:
            self.value = default
        self.name = name


class ParMeta(type):
    def __new__(cls, name, bases, clsdict):
        clsdict['_parms'] = [key for key, val in clsdict.items()
            if isinstance(val, Parm)]
        print(clsdict['_parms'])
        clsobj = super().__new__(cls, name, bases, clsdict)
        return clsobj

class Proto(metaclass=ParMeta):
    x = Parm(1.)
    y = Parm(2.5)
    def __init__(self):
        pass

    def __setattr__(self, name, value):
        if name in self._parms:
            self.__class__.__dict__[name].value = value
        else:
            self.__dict__[name] = value


