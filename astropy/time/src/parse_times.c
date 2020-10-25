#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include <stdio.h>
#include <string.h>

#define MODULE_DOCSTRING \
    "Fast time parsers.\n\n" \
    "This module allows one to create gufuncs that vectorize the parsing of\n" \
    "standard time strings."
#define CREATE_PARSER_DOCSTRING \
    "create_parser()\n\n" \
    "Create a gufunc that parsers according to the passed in parameters.\n\n" \
    "Parameters\n" \
    "----------\n" \
    "pars : ~numpy.ndarray\n" \
    "    Should be structured array with delim, start, stop, break_allowed for each\n" \
    "    of year, month, day, hour, minute, integer second, fractional second.\n\n" \
    "Returns\n" \
    "-------\n" \
    "parser : `~numpy.ufunc`\n" \
    "    Will parse byte strings passed in as `~numpy.uint8`."


// ASCII codes for '0' and '9'
const char char_zero = 48;
const char char_nine = 57;

struct pars_struct_t {
    char delim;
    int start;
    int stop;
    npy_bool break_allowed;
};

// Distutils on Windows automatically exports ``PyInit__parse_times``,
// create dummy to prevent linker complaining about missing symbol.
// Based on convolution/src/convolve.c.
#if defined(_MSC_VER)
void PyInit__parse_times(void)
{
    return;
}
#endif

int parse_int_from_char_array(char *chars, int str_len,
                              char delim, int idx0, int idx1,
                              int *val)
// Parse integer from positions idx0:idx1 (inclusive) within chars, optionally
// starting with a delimiter.
//
// Example: "2020-01-24"
//                  ^^^
//           0123456789
//
// int day, status;
// status = parse_int_from_char_array("2020-01-24", 10, '-', 7, 9, &day);
//
// Inputs:
//  char *chars: time string
//  int str_len: length of *chars string
//  char delim: optional character at position idx0 when delim > 0
//  int idx0: start index for parsing integer
//  int idx1: stop index (inclusive) for parsing integer
//
// Output:
//  int *val: output value
//
// Returns:
//  int status:
//    0: OK
//    1: String ends at the beginning of requested value
//    2: String ends in the middle of requested value
//    3: Required delimiter character not found
//    4: Non-digit found where digit (0-9) required
{
    int mult = 1;
    char digit;
    char ch;
    int ii;

    // Check if string ends (has 0x00) before str_len. Require that this segment
    // of the string is entirely contained in the string (idx1 < str_len),
    // remembering that idx1 is inclusive and counts from 0.
    if (idx1 < str_len) {
        for (ii = idx0; ii <= idx1; ii++) {
            if (chars[ii] == 0) {
                str_len = ii;
                break;
            }
        }
    }
    // String ends before the beginning of requested value,
    // e.g. "2000-01" (str_len=7) for day (idx0=7). This is OK in some
    // cases, e.g. before hour (2000-01-01).
    if (idx0 >= str_len) {
        return 1;
    }

    // String ends in the middle of requested value. This implies a badly
    // formatted time.
    if (idx1 >= str_len) {
        return 2;
    }

    // Look for optional delimiter character, e.g. ':' before minute. If delim == 0
    // then no character is required.
    if (delim > 0) {
        // Required start character not found.
        if (chars[idx0] != delim) {
            return 3;
        }
        idx0 += 1;
    }

    // Build up the value using reversed digits
    *val = 0;
    for (ii = idx1; ii >= idx0; ii--)
    {
        ch = chars[ii];
        if (ch < char_zero || ch > char_nine) {
            // Not a digit, implying badly formatted time.
            return 4;
        }
        digit = ch - char_zero;
        *val += digit * mult;
        mult *= 10;
    }

    return 0;
}

int parse_frac_from_char_array(char *chars, int str_len, char delim, int idx0,
                               double *val)
// Parse trailing fraction starting from position idx0 in chars.
//
// Example: "2020-01-24T12:13:14.5556"
//                              ^^^^^
//           012345678901234567890123
//
// int status;
// float frac;
// status = parse_frac_from_char_array("2020-01-24T12:13:14.5556", 24, '.', 19, &frac);
//
// Inputs:
//  char *chars: time string
//  int str_len: length of *chars string
//  char delim: optional character at position idx0 when delim > 0
//  int idx0: start index for parsing integer
//
// Output:
//  double *val: output value
//
// Returns:
//  int status:
//    0: OK
//    1: String ends at the beginning of requested value
//    3: Required delimiter character not found
//    4: Non-digit found where digit (0-9) required
{
    double mult = 0.1;
    char digit;
    char ch;
    int ii;

    *val = 0.0;

    // String ends at exactly before the beginning of requested fraction.
    // e.g. "2000-01-01 12:13:14". Fraction value is zero.
    if (idx0 == str_len) {
        return 1;
    }

    // Look for optional delimiter character, e.g. '.' before fraction. If delim == 0
    // then no character is required. This can happen for unusual formats like
    // Chandra GRETA time yyyyddd.hhmmssfff.
    if (delim > 0) {
        // Required start character not found.
        if (chars[idx0] != delim) {
            return 3;
        }
        idx0 += 1;
    }

    for (ii = idx0; ii < str_len; ii++)
    {
        ch = chars[ii];
        if (ch < char_zero || ch > char_nine) {
            // Not a digit, implying badly formatted time.
            return 4;
        }
        digit = ch - char_zero;
        *val += digit * mult;
        mult /= 10.0;
    }
    return 0;
}

static inline int is_leap_year (int year)
// Determine if year is a leap year.
// Inspired by from https://stackoverflow.com/questions/17634282
{
  return ((year & 3) == 0)
          && ((year % 100 != 0)
              || (((year / 100) & 3) == 0));
}

int convert_day_of_year_to_month_day(int year, int day_of_year, int *month, int *day_of_month)
// Convert year and day_of_year into month, day_of_month
// Inspired by from https://stackoverflow.com/questions/17634282, determine
{
    int leap_year = is_leap_year(year) ? 1 : 0;
    int days_in_year = leap_year ? 366 : 365;
    const unsigned short int _mon_yday_normal[13] =
        { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
    const unsigned short int _mon_yday_leap[13] =
        { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };
    const unsigned short int *mon_yday = leap_year ? _mon_yday_leap :_mon_yday_normal;
    int mon;

    if (day_of_year < 1 || day_of_year > days_in_year) {
        // Error in day_of_year
        return 5;
    }

    for (mon = 1; mon <= 12; mon++) {
        if (day_of_year <= mon_yday[mon]) {
            *month = mon;
            *day_of_month = day_of_year - mon_yday[mon - 1];
            break;
        }
    }

    return 0;
}


static void
parser_loop(char **args, npy_intp *dimensions, npy_intp *steps, void *data)
{
    npy_intp n = dimensions[0];
    npy_intp max_str_len = dimensions[1];

    struct pars_struct_t *pars = (struct pars_struct_t *)PyArray_DATA((PyArrayObject *)data);
    npy_bool has_day_of_year = pars[1].start < 0;

    npy_intp ii, i, str_len, status;

    static char *msgs[5] = {
        "time string ends at beginning of component where break is not allowed",
        "time string ends in middle of component",
        "required delimiter character not found",
        "non-digit found where digit (0-9) required",
        "bad day of year (1 <= doy <= 365 or 366 for leap year"};

    for (ii = 0; ii < n; ii++)
    {
        char *time = args[0];
        char **out = args+1;

        // Check for null termination before max_str_len. If called using a contiguous
        // numpy 2-d array of chars there may or may not be null terminations.
        str_len = max_str_len;
        for (i = 0; i < max_str_len; i++) {
            if (time[i] == 0) {
                str_len = i;
                break;
            }
        }
        // Parse "2000-01-12 13:14:15.678"
        //        01234567890123456789012
        // pars is a 7-element struct with for each of year, month, day,
        // hour, minute, second_int, and second_frac containing the delimiter
        // (if any) that starts it, start and stop index of the relevant part
        // of the string, and whether a break is allowed before this item.

        // Here, first get the year, month, day, hour and minute parts,
        // store them directly into the relevant integer output arrays.
        for (i=0; i < 5; i++) {
            if (i == 1 && has_day_of_year) continue;
            status = parse_int_from_char_array(
                time, str_len, pars[i].delim, pars[i].start, pars[i].stop,
                (int *)out[i]);
            if (status) break;
        }
        if (status) {
            // Break or error.
            if (!(status == 1 && pars[i].break_allowed)) goto error;
            // If we had an allowed break, i.e., a short but correct string,
            // set the output arrays to defaults for the remaining elements.
            while (i < 3) {
                *(int *)out[i++] = 1;  // month, day.
            }
            while (i < 5) {
                *(int *)out[i++] = 0;  // hour, minute.
            }
            *(double *)out[5] = 0.;    // second.
        }
        else {
            // No break so far; parse integer and fractional parts of seconds.
            int second_int = 0;
            double second_frac = 0.;

            status = parse_int_from_char_array(
                time, str_len, pars[i].delim, pars[i].start, pars[i].stop,
                &second_int);
            if (!status) {
                i++;
                status = parse_frac_from_char_array(
                    time, str_len, pars[i].delim, pars[i].start, &second_frac);
            }
            if (status && !(status == 1 && pars[i].break_allowed)) goto error;
            // Fill output for seconds with the sum of integer and fraction parts.
            *(double *)out[5] = (double)second_int + second_frac;
        }

        if (has_day_of_year) {
            // day contains day of year at this point, while month is unset.
            // convert it to month and day of month (with year passed on as
            // well to allow checks for leap years).
            status = convert_day_of_year_to_month_day(
                *(int *)out[0], *(int *)out[2], (int *)out[1], (int *)out[2]);
            if (status) { goto error; }
        }

        for (i=0; i<7; i++) {
            args[i] += steps[i];
        }

    }
    return;

  error:
    PyErr_Format(PyExc_ValueError,
                 "fast C time string parser failed: %s", msgs[status-1]);
    return;
}


/* Create a gufunc parser */

static PyArray_Descr *dt_pars = NULL;   /* Set in PyInit_ufunc */
static char gufunc_types[] = {
    NPY_UINT8, NPY_INT32, NPY_INT32, NPY_INT32,
    NPY_INT32, NPY_INT32, NPY_DOUBLE};
static PyUFuncGenericFunction parser_loops[] = {parser_loop};

static PyObject *
create_parser(PyObject *NPY_UNUSED(dummy), PyObject *args, PyObject *kwds)
{
    /* Input arguments */
    char *kw_list[] = {"pars", "name", "doc", NULL};
    PyObject *pars;
    char *name=NULL, *doc=NULL;
    /* Output */
    PyUFuncObject *gufunc=NULL;

    PyArrayObject *pars_array;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|ss", kw_list,
                                     &pars, &name, &doc)) {
        return NULL;
    }
    if (name == NULL) {
        name = "fast_parser";
    }
    Py_INCREF(dt_pars);
    pars_array = (PyArrayObject *)PyArray_FromAny(pars, dt_pars, 1, 1,
                     (NPY_ARRAY_CARRAY | NPY_ARRAY_ENSURECOPY), NULL);
    if (pars_array == NULL) {
        return NULL;
    }
    if (PyArray_SIZE(pars_array) != 7) {
        PyErr_SetString(PyExc_ValueError,
                        "Parameter array must have 7 entries"
                        "(year, month, day, hour, minute, integer second, fraction)");
    }

    gufunc = (PyUFuncObject *)PyUFunc_FromFuncAndDataAndSignature(
        parser_loops, NULL, gufunc_types, 1,
        1, 6, PyUFunc_None, name, doc, 0, "(n)->(),(),(),(),(),()");
    if (gufunc == NULL) {
        goto fail;
    }
    /*
     * We need to keep array around, as this has the required information, but
     * it should be deallocated when the ufunc is deleted. Use ->obj for this
     * (also used in frompyfunc).
     */
    gufunc->obj = pars_array;
    gufunc->data = &gufunc->obj;
    return (PyObject *)gufunc;

  fail:
    Py_XDECREF(pars_array);
    Py_XDECREF(gufunc);
    return NULL;
}


static PyMethodDef parse_times_methods[] = {
    {"create_parser", (PyCFunction)create_parser,
         METH_VARARGS | METH_KEYWORDS, CREATE_PARSER_DOCSTRING},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "parse_times",
        MODULE_DOCSTRING,
        -1,
        parse_times_methods,
        NULL,
        NULL,
        NULL,
        NULL
};

/* Initialization function for the module */
PyMODINIT_FUNC PyInit__parse_times(void) {
    PyObject *m;
    PyObject *d;
    PyObject *dtype_def;

    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }
    import_array();
    import_ufunc();

    d = PyModule_GetDict(m);

    /* parameter and output dtypes */
    dtype_def = Py_BuildValue(
        "[(s, s), (s, s), (s, s), (s, s)]",
        "delim", "S1",
        "start", "i4",
        "stop", "i4",
        "break_allowed", "?");
    PyArray_DescrAlignConverter(dtype_def, &dt_pars);
    Py_DECREF(dtype_def);

    if (dt_pars == NULL) {
        goto fail;
    }
    PyDict_SetItemString(d, "dt_pars", (PyObject *)dt_pars);

    goto decref;

  fail:
    Py_XDECREF(m);
    m = NULL;

  decref:
    Py_XDECREF(dt_pars);
    return m;
}
