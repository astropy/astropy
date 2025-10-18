# Licensed under a 3-clause BSD style license - see LICENSE.rst


from astropy.cosmology._src.parameter import Parameter, all_parameters


def test_all_parameters():
    """Test :func:`astropy.cosmology._src.utils.all_parameters`."""

    class ClassA:
        a = 1
        b = 2

    got = all_parameters(ClassA)
    assert got == {}

    class ClassB(ClassA):
        H0: Parameter = Parameter(
            doc="Hubble constant at z=0.",
            unit="km/(s Mpc)",
            fvalidate="scalar",
        )

    got = all_parameters(ClassB)
    assert got.keys() == {"H0"}
    assert all(isinstance(p, Parameter) for p in got.values())
