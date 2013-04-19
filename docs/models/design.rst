.. _models-design:

*******************
Models Design Goals
*******************

The `~astropy.models.models` and `~astropy.models.fitting` modules described
here are designed to work as peers with each other. The goal is to be able to
add models without explicit reference to fitting algorithms and likewise, add 
different fitting algorithms without changing the existing models. The mechanism 
that allows this is the special `~astropy.models.parameters` module that both models 
and fitters use to interact with each other. Nevertheless, most users won't need to 
interact with this module unless they wish to add new models or 
fitters (the term used hereafter for specific fitting algorithms) to 
the existing suites of models and fitters.

Furthermore, the models are designed to be combined in many ways. It
is possible, for example, to combine models serially
`~astropy.models.models.SCompositeModel`, so that the output values of one model are
used as input values to another. It is also possible to form a new model by
combining models in parallel (each model is evaluated separately with the
original input and the deltas are summed), `~astropy.models.models.PCompositeModel`.
Since models may have multiple input values, machinery is provided that allows
assigning outputs from one model into the appropriate input of another in a
flexible way, `~astropy.models.models.LabeledInput`. Finally, it is permitted
to combine any number of models using all of these mechanisms simultaneously.
A composite model can be used to make further composite models.

In the future this will support a model language which will allow using models
in algebraic operations like

.. math:: model = (model_1 + model_2) * model_3
