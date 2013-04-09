.. _parameters:

**********
Parameters
**********

Parameters are used in three different contexts within this package: 
communicating with fitters, model evaluation and getting values to/from users. 

Models maintain a list of parameter names, **`parnames`**. Single parameters are list-like 
objects, instances of **`parameters._Parameter`**. Simple mathematical operations can be 
performed with them. The preferred way for users to interact with models is through 
individual parameters.

The goal of this package is, when possible, to allow simultaneous model evaluation 
and fitting with multiple parameter sets. Because of this, all models have a **`psets`**
attribute, an array of shape `(len(parnames), paramdim)`, where **`paramdim`** is the number of 
parameter sets. Typically the array is of type float but can become an object array in 
some cases. **`Psets`** is used for model evaluation.

In addition, all models maintain an attribute,
**`parameters`**, an instance of **`parameters.Parameters`**. This is a flat list of 
parameter values which fitters update. It serves as a communication tool between fitters 
and models.

Individual parameters, **`psets`** and the flat list of parameter values are kept in sync. 
Single parameters are updated through properties. An update to a single parameter 
trigers an update to **`psets`** and **`parameters`**. Single parameters are updated 
after a change to **`parameters`**. **`Psets`** are always constructed on demand from single 
parameters and cannot be updated directly.

Parameters Examples
-------------------

- Polynomial models are created by default with all coefficients set to 0.

>>> p1=models.Poly1DModel(degree=4)
>>> p1.parnames
['c0', 'c1', 'c2', 'c3', 'c4']
>>> p1.parameters
[0.0, 0.0, 0.0, 0.0, 0.0]

- Coefficients can be set using the **`parameters`** attribute

>>> p1.parameters = [0, 1, 2, 3, 4]
>>> p1.parameters
[0.0, 1.0, 2.0, 3.0, 4.0]

- It is possible to set the coefficients passing a dictionary to the `pars` keyword

>>> ch2 = models.ICheb2DModel(xdeg=2,ydeg=3, paramdim=2)
>>> coeff = {}
>>> for i, j in zip(ch2.parnames, range(len(ch2.parnames))):
        coeff[i] = [j, j+10]
>>> ch2 = models.ICheb2DModel(xdeg=2, ydeg=3, **coeff)
>>> ch2.psets
array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11],
       [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]])


- or directly, using keyword arguments

>>> ch2 = models.ICheb2DModel(xdeg=2, ydeg=3, 'c0_0'=[0, 10], 'c0_1'=[3, 13],
        ...                                       'c0_2'=[6, 16], 'c0_3'=[9, 19],
        ...                                       'c1_0'=[1, 11], 'c1_1'=[4, 14],
        ...                                       'c1_2'=[7, 17], 'c1_3'=[10, 20,],
        ...                                       'c2_0'=[2, 12], 'c2_1'=[5, 15],
        ...                                       'c2_2'=[8, 18], 'c2_3'=[11, 21])


- It is possible to change a single parameter

>>> ch2.c0_0
[0, 100]
>>> ch2.c0_0[0] = -34.2
>>> ch2.c0_0
[-34.2, 10]

- The number of parameter sets is stored in an attribute `paramdim`

>>> ch2.paramdim
2


    