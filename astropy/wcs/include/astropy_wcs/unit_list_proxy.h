/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#ifndef __UNIT_LIST_PROXY_H__
#define __UNIT_LIST_PROXY_H__

#include "pyutil.h"

/***************************************************************************
 * List-of-units proxy object
 *
 * A Python object that looks like a list of units, but is back by a C
 *   char * list[];
 ***************************************************************************/

/*@null@*/ PyObject *
PyUnitListProxy_New(
    PyObject* owner,
    Py_ssize_t size,
    char (*array)[72]
    );

int
_setup_unit_list_proxy_type(
    PyObject* m);

static INLINE PyObject*
get_unit_list(
    /*@unused@*/ const char* propname,
    char (*array)[72],
    Py_ssize_t len,
    PyObject* owner) {

  return PyUnitListProxy_New(owner, len, array);
}

int
set_unit_list(
    PyObject *owner,
    const char* propname,
    PyObject* value,
    Py_ssize_t len,
    char (*dest)[72]);

#endif /* __UNIT_LIST_PROXY_H__ */
