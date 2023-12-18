# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .codegen import make_function_with_signature as make_function_with_signature
from .decorators import (
    classproperty as classproperty,
    deprecated as deprecated,
    deprecated_attribute as deprecated_attribute,
    deprecated_renamed_argument as deprecated_renamed_argument,
    format_doc as format_doc,
    lazyproperty as lazyproperty,
    sharedmethod as sharedmethod,
)
from .introspection import (
    find_current_module as find_current_module,
    isinstancemethod as isinstancemethod,
    minversion as minversion,
    resolve_name as resolve_name,
)
from .misc import (
    JsonCustomEncoder as JsonCustomEncoder,
    NumpyRNGContext as NumpyRNGContext,
    dtype_bytes_or_chars as dtype_bytes_or_chars,
    find_api_page as find_api_page,
    format_exception as format_exception,
    indent as indent,
    is_path_hidden as is_path_hidden,
    isiterable as isiterable,
    silence as silence,
    walk_skip_hidden as walk_skip_hidden,
)
from .shapes import (
    IncompatibleShapeError as IncompatibleShapeError,
    NDArrayShapeMethods as NDArrayShapeMethods,
    ShapedLikeNDArray as ShapedLikeNDArray,
    check_broadcast as check_broadcast,
    simplify_basic_index as simplify_basic_index,
    unbroadcast as unbroadcast,
)
from . import (
    argparse as argparse,
    codegen as codegen,
    collections as collections,
    data_info as data_info,
    diff as diff,
    exceptions as exceptions,
    introspection as introspection,
    misc as misc,
    parsing as parsing,
    setup_package as setup_package,
    shapes as shapes,
    state as state,
    _compiler as _compiler,
    console as console,
    data as data,
    decorators as decorators,
    compat as compat,
    iers as iers,
    masked as masked,
    metadata as metadata,
    tests as tests,
    xml as xml,
)
