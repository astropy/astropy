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
    "    Should have delims, start, stop, break_allowed for each of\n" \
    "    year, month, day, hour, minute, integer second, fractional second.\n\n" \
    "Returns\n" \
    "-------\n" \
    "parser : `~numpy.ufunc`\n" \
    "    Will parser bytes or unicode."


// ASCII codes for '0' and '9'
const char char_zero = 48;
const char char_nine = 57;

struct time_struct_t {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second_int;
    double second_frac;
};

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
// status = parse_int_from_char_array("2020-01-24", &day, 10, '-', 7, 9);
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
// status = parse_frac_from_char_array("2020-01-24T12:13:14.5556", &frac, 24, '.', 19);
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
    char *time = args[0];
    char *tm_ptr = args[1];
    npy_intp i_time=steps[0];
    npy_intp i_tm=steps[1];
    struct pars_struct_t *pars;

    npy_intp ii, i, str_len, status;
    struct time_struct_t *tm;
    npy_bool has_day_of_year;

    static char *msgs[5] = {
        "time string ends at beginning of component where break is not allowed",
        "time string ends in middle of component",
        "required delimiter character not found",
        "non-digit found where digit (0-9) required",
        "bad day of year (1 <= doy <= 365 or 366 for leap year"};

    for (ii = 0; ii < n; ii++, time+=i_time, tm_ptr+=i_tm)
    {
        // Initialize default values
        tm = (struct time_struc_t *)tm_ptr;
        tm->month = 1;
        tm->day = 1;
        tm->hour = 0;
        tm->minute = 0;
        tm->second_int = 0;
        tm->second_frac = 0.0;

        // Parse "2000-01-12 13:14:15.678"
        //        01234567890123456789012

        // Check for null termination before max_str_len. If called using a contiguous
        // numpy 2-d array of chars there may or may not be null terminations.
        str_len = max_str_len;
        for (i = 0; i < max_str_len; i++) {
            if (time[i] == 0) {
                str_len = i;
                break;
            }
        }

        // Get each time component year, month, day, hour, minute, isec, frac
        pars = (struct pars_struct_t *)data;
        status = parse_int_from_char_array(time, str_len, pars->delim, pars->start, pars->stop, &tm->year);
        if (status) {
            if (status == 1 && pars->break_allowed) { continue; }
            else { goto error; }
        }

        pars++;
        has_day_of_year = pars->start < 0;
        // Optionally parse month
        if (!has_day_of_year) {
            status = parse_int_from_char_array(time, str_len, pars->delim, pars->start, pars->stop, &tm->month);
            if (status) {
                if (status == 1 && pars->break_allowed) { continue; }
                else { goto error; }
            }
        }

        pars++;
        // This might be day-of-month or day-of-year
        status = parse_int_from_char_array(time, str_len, pars->delim, pars->start, pars->stop, &tm->day);
        if (status) {
            if (status == 1 && pars->break_allowed) { continue; }
            else { goto error; }
        }

        if (has_day_of_year) {
            // day contains day of year at this point, but convert it to day of month
            status = convert_day_of_year_to_month_day(tm->year, tm->day, &tm->month, &tm->day);
            if (status) { goto error; }
        }

        pars++;
        status = parse_int_from_char_array(time, str_len, pars->delim, pars->start, pars->stop, &tm->hour);
        if (status) {
            if (status == 1 && pars->break_allowed) { continue; }
            else { goto error; }
        }

        pars++;
        status = parse_int_from_char_array(time, str_len, pars->delim, pars->start, pars->stop, &tm->minute);
        if (status) {
            if (status == 1 && pars->break_allowed) { continue; }
            else { goto error; }
        }

        pars++;
        status = parse_int_from_char_array(time, str_len, pars->delim, pars->start, pars->stop, &tm->second_int);
        if (status) {
            if (status == 1 && pars->break_allowed) { continue; }
            else { goto error; }
        }

        pars++;
        status = parse_frac_from_char_array(time, str_len, pars->delim, pars->start, &tm->second_frac);
        if (status) {
            if (status != 1 || !pars->break_allowed) { goto error; }
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
static PyArray_Descr *gufunc_dtypes[2];


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
    int status;

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
        NULL, NULL, NULL, 0, 1, 1, PyUFunc_None, name, doc, 0, "(n)->()");
    if (gufunc == NULL) {
        goto fail;
    }
    status = PyUFunc_RegisterLoopForDescr(
        gufunc, gufunc_dtypes[0], parser_loop, gufunc_dtypes, PyArray_DATA(pars_array));
    if (status != 0) {
        goto fail;
    }
    /*
     * We need to keep array around, as this has the required information, but
     * it should be deallocated when the ufunc is deleted. Use ->obj for this
     * (also used in frompyfunc).
     */
    gufunc->obj = pars_array;
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
    PyArray_Descr *dt_u1 = NULL, *dt_ymdhmsf = NULL;

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

    dtype_def = Py_BuildValue("[(s, s)]", "byte", "u1");
    PyArray_DescrAlignConverter(dtype_def, &dt_u1);
    Py_DECREF(dtype_def);

    dtype_def = Py_BuildValue(
        "[(s, s), (s, s), (s, s), (s, s), (s, s), (s, s), (s, s)]",
        "year", "i4",
        "month", "i4",
        "day", "i4",
        "hour", "i4",
        "minute", "i4",
        "second_int", "i4",
        "second_frac", "f8");
    PyArray_DescrAlignConverter(dtype_def, &dt_ymdhmsf);
    Py_DECREF(dtype_def);
    if (dt_pars == NULL || dt_u1 == NULL || dt_ymdhmsf == NULL) {
        goto fail;
    }
    PyDict_SetItemString(d, "dt_pars", (PyObject *)dt_pars);
    PyDict_SetItemString(d, "dt_u1", (PyObject *)dt_u1);
    PyDict_SetItemString(d, "dt_ymdhmsf", (PyObject *)dt_ymdhmsf);

    gufunc_dtypes[0] = dt_u1;
    gufunc_dtypes[1] = dt_ymdhmsf;

    goto decref;

  fail:
    Py_XDECREF(m);
    m = NULL;

  decref:
    Py_XDECREF(dt_pars);
    Py_XDECREF(dt_u1);
    Py_XDECREF(dt_ymdhmsf);
    return m;
}
