/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __STR_LIST_PROXY_H__
#define __STR_LIST_PROXY_H__

#include "pyutil.h"

/***************************************************************************
 * List-of-strings proxy object
 *
 * A Python object that looks like a list of strings, but is back by a C
 *   char * list[];
 ***************************************************************************/

typedef int (*str_verify_fn)(const char *);

/*@null@*/ PyObject *
PyStrListProxy_New(
    PyObject* owner,
    Py_ssize_t size,
    Py_ssize_t maxsize,
    char (*array)[72]
    );

/*@null@*/ PyObject*
str_list_proxy_repr(
    char (*array)[72],
    Py_ssize_t size,
    Py_ssize_t maxsize);

int
_setup_str_list_proxy_type(
    PyObject* m);

#endif /* __STR_LIST_PROXY_H__ */
