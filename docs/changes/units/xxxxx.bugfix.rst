``Quantity.__array_ufunc__()`` now returns ``NotImplemented`` instead of raising
``ValueError`` when it cannot convert inputs, allowing duck array types to handle
operations via their reflected operators. This enables better interoperability
with custom array types that implement ``__array_ufunc__`` and have unit handling.
