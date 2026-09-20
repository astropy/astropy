# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Helpers for the LALR(1) parsers (built with `lark`) used to parse unit
and angle strings.
"""

from lark import Lark, Token, Transformer

__all__ = []


def token_values(f, _data, children, _meta):
    """Call a transformer method with the values of the child tokens.

    Meant to be used as ``@lark.v_args(wrapper=token_values)`` on a
    `lark.Transformer` subclass: rule methods then receive the *values* of
    terminal tokens (as set by the ``t_<TERMINAL>`` methods, see
    `make_parser`) rather than the `lark.Token` instances, and the results
    of sub-rules, as positional arguments.
    """
    return f(
        *[child.value if isinstance(child, Token) else child for child in children]
    )


def make_parser(grammar: str, transformer: Transformer, start: str) -> Lark:
    """Create an LALR(1) parser from a grammar and a transformer.

    The parser applies the transformer callbacks while parsing, so that
    ``parser.parse(text)`` directly returns the transformed result.
    Methods of the transformer called ``t_<TERMINAL>`` are applied to tokens
    of that terminal as soon as they are produced by the lexer (like the
    rules of a ``lex`` lexer), before any grammar rule is applied to them.
    They should set the ``value`` of the token and return it.

    The returned parser is stateless and can safely be shared between threads.
    """
    lexer_callbacks = {
        name.removeprefix("t_"): getattr(transformer, name)
        for name in dir(transformer)
        if name.startswith("t_")
    }
    return Lark(
        grammar,
        parser="lalr",
        lexer="basic",
        transformer=transformer,
        lexer_callbacks=lexer_callbacks,
        start=start,
    )
