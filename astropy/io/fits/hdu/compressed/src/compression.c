#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <fits_hcompress.h>
#include <fits_hdecompress.h>
#include <pliocomp.h>
#include <quantize.h>
#include <unquantize.h>
#include <ricecomp.h>


// Compatibility code because we pick up fitsio2.h from cextern. Can
// remove once we remove cextern
#ifdef _REENTRANT
pthread_mutex_t Fitsio_Lock;
int Fitsio_Pthread_Status = 0;
#endif

/* Define docstrings */
static char module_docstring[] = "Core compression/decompression functions wrapped from cfitsio.";
static char compress_plio_1_c_docstring[] = "Compress data using PLIO_1";
static char decompress_plio_1_c_docstring[] = "Decompress data using PLIO_1";
static char compress_rice_1_c_docstring[] = "Compress data using RICE_1";
static char decompress_rice_1_c_docstring[] = "Decompress data using RICE_1";
static char compress_hcompress_1_c_docstring[] =
    "Compress data using HCOMPRESS_1";
static char decompress_hcompress_1_c_docstring[] =
    "Decompress data using HCOMPRESS_1";
static char quantize_float_c_docstring[] = "Quantize float data";
static char quantize_double_c_docstring[] = "Quantize float data";
static char unquantize_float_c_docstring[] = "Unquantize data to float";
static char unquantize_double_c_docstring[] = "Unquantize data to double";

/* Declare the C functions here. */
static PyObject *compress_plio_1_c(PyObject *self, PyObject *args);
static PyObject *decompress_plio_1_c(PyObject *self, PyObject *args);
static PyObject *compress_rice_1_c(PyObject *self, PyObject *args);
static PyObject *decompress_rice_1_c(PyObject *self, PyObject *args);
static PyObject *compress_hcompress_1_c(PyObject *self, PyObject *args);
static PyObject *decompress_hcompress_1_c(PyObject *self, PyObject *args);
static PyObject *quantize_float_c(PyObject *self, PyObject *args);
static PyObject *quantize_double_c(PyObject *self, PyObject *args);
static PyObject *unquantize_float_c(PyObject *self, PyObject *args);
static PyObject *unquantize_double_c(PyObject *self, PyObject *args);
static PyObject *CfitsioException = NULL;

/* Define the methods that will be available on the module. */
static PyMethodDef module_methods[] = {
    {"compress_plio_1_c", compress_plio_1_c, METH_VARARGS, compress_plio_1_c_docstring},
    {"decompress_plio_1_c", decompress_plio_1_c, METH_VARARGS, decompress_plio_1_c_docstring},
    {"compress_rice_1_c", compress_rice_1_c, METH_VARARGS, compress_rice_1_c_docstring},
    {"decompress_rice_1_c", decompress_rice_1_c, METH_VARARGS, decompress_rice_1_c_docstring},
    {"compress_hcompress_1_c", compress_hcompress_1_c, METH_VARARGS, compress_hcompress_1_c_docstring},
    {"decompress_hcompress_1_c", decompress_hcompress_1_c, METH_VARARGS, decompress_hcompress_1_c_docstring},
    {"quantize_float_c", quantize_float_c, METH_VARARGS, quantize_float_c_docstring},
    {"quantize_double_c", quantize_double_c, METH_VARARGS, quantize_double_c_docstring},
    {"unquantize_float_c", unquantize_float_c, METH_VARARGS, unquantize_float_c_docstring},
    {"unquantize_double_c", unquantize_double_c, METH_VARARGS, unquantize_double_c_docstring},
    {NULL, NULL, 0, NULL}
};

/* This is the function that is called on import. */
static PyModuleDef compression = {
    PyModuleDef_HEAD_INIT,
    "_compression",
    module_docstring,
    -1,
    module_methods,
};

PyMODINIT_FUNC
PyInit__compression(void)
{
    PyObject* m;
    m = PyModule_Create(&compression);

    /* Initialize new exception object */
    CfitsioException = PyErr_NewException("_compression.CfitsioException", NULL, NULL);

    /* Add exception object to your module */
    PyModule_AddObject(m, "CfitsioException", CfitsioException);

    return m;
};

// Some of the cfitsio compression functions use this function to put an error
// message on the stack.  We provide an implementation which sets the Python
// error state with our custom exception type, and the message provided by the
// cfitsio call.  In our wrapper functions we can then check if the Python error
// state is set and then return NULL to raise the error.
void ffpmsg(const char *err_message) {
    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();
    PyErr_SetString(CfitsioException, err_message);
    PyGILState_Release(gstate);
}

/* PLIO/IRAF compression */

static PyObject *compress_plio_1_c(PyObject *self, PyObject *args) {

  const char *str;
  char *buf;
  Py_ssize_t count;
  PyObject *result;

  int maxelem;
  int tilesize;
  short *compressed_values;
  int compressed_length;
  int *decompressed_values;

  if (!PyArg_ParseTuple(args, "y#i", &str, &count, &tilesize)) {
    return NULL;
  }

  decompressed_values = (int *)str;

  for (int ii = 0; ii < tilesize; ii++)  {
    if (decompressed_values[ii] < 0 || decompressed_values[ii] > 16777215)
    {
      /* plio algorithm only supports positive 24 bit ints */
      PyErr_SetString(PyExc_ValueError,
                      "data out of range for PLIO compression (0 - 2**24)");
      return (PyObject *)NULL;
    }
  }

  // For PLIO imcomp_calc_max_elem in cfitsio does this to calculate max memory:
  maxelem = tilesize;
  // However, when compressing small numbers of random integers you can end up
  // using more memory for the compressed bytes.  In the worst case scenario we
  // tested, compressing a single 4 byte integer will compress to 16 bytes.  We
  // therefore allocate a buffer 4 ints larger than we need here to give that
  // margin of error.
  compressed_values = (short *)calloc(maxelem + 4, sizeof(int));

  decompressed_values = (int *)str;

  compressed_length = pl_p2li(decompressed_values, 1, compressed_values, tilesize);

  if (PyErr_Occurred() != NULL) {
    // If an error condition inside the cfitsio function, the call inside
    // cfitsio should have called the ffpmsg function which sets the Python
    // exception, so we just return here to raise an error.
    return (PyObject *)NULL;
  }

  buf = (char *)compressed_values;

  result = Py_BuildValue("y#", buf, compressed_length * 2);
  free(buf);
  return result;
}

static PyObject *decompress_plio_1_c(PyObject *self, PyObject *args) {

  const char *str;
  char *buf;
  Py_ssize_t count;
  PyObject *result;

  int tilesize;

  short *compressed_values;
  int *decompressed_values;

  if (!PyArg_ParseTuple(args, "y#i", &str, &count, &tilesize)) {
    return NULL;
  }

  compressed_values = (short *)str;

  decompressed_values = (int *)calloc(tilesize, sizeof(int));

  pl_l2pi(compressed_values, 1, decompressed_values, tilesize);

  if (PyErr_Occurred() != NULL) {
    // If an error condition inside the cfitsio function, the call inside
    // cfitsio should have called the ffpmsg function which sets the Python
    // exception, so we just return here to raise an error.
    return (PyObject *)NULL;
  }

  buf = (char *)decompressed_values;

  result = Py_BuildValue("y#", buf, tilesize * sizeof(int));
  free(buf);
  return result;
}

/* RICE compression */

static PyObject *compress_rice_1_c(PyObject *self, PyObject *args) {

  const char *str;
  Py_ssize_t count;
  PyObject *result;

  int blocksize, bytepix;

  int maxelem;
  unsigned char *compressed_values;
  int compressed_length;
  signed char *decompressed_values_byte;
  short *decompressed_values_short;
  int *decompressed_values_int;

  if (!PyArg_ParseTuple(args, "y#ii", &str, &count, &blocksize, &bytepix)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS

  // maxelem adapted from cfitsio's imcomp_calc_max_elem function
  maxelem = count + count / bytepix / blocksize + 2 + 4;

  compressed_values = (unsigned char *)malloc(maxelem);

  if (bytepix == 1) {
    decompressed_values_byte = (signed char *)str;
    compressed_length = fits_rcomp_byte(decompressed_values_byte, (int)count, compressed_values, count * 16, blocksize);
  } else if (bytepix == 2) {
    decompressed_values_short = (short *)str;
    compressed_length = fits_rcomp_short(decompressed_values_short, (int)count / 2, compressed_values, count * 16, blocksize);
  } else {
    decompressed_values_int = (int *)str;
    compressed_length = fits_rcomp(decompressed_values_int, (int)count / 4, compressed_values, count * 16, blocksize);
  }

  Py_END_ALLOW_THREADS

  if (PyErr_Occurred() != NULL) {
    // If an error condition inside the cfitsio function, the call inside
    // cfitsio should have called the ffpmsg function which sets the Python
    // exception, so we just return here to raise an error.
    return (PyObject *)NULL;
  }

  result = Py_BuildValue("y#", compressed_values, compressed_length);
  free(compressed_values);
  return result;
}

static PyObject *decompress_rice_1_c(PyObject *self, PyObject *args) {

  const char *str;
  char *dbytes;
  Py_ssize_t count;
  PyObject *result;

  int blocksize, bytepix, tilesize;

  unsigned char *compressed_values;
  unsigned char *decompressed_values_byte;
  unsigned short *decompressed_values_short;
  unsigned int *decompressed_values_int;

  if (!PyArg_ParseTuple(args, "y#iii", &str, &count, &blocksize, &bytepix, &tilesize)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS

  compressed_values = (unsigned char *)str;

  if (bytepix == 1) {
    decompressed_values_byte = (unsigned char *)malloc(tilesize);
    fits_rdecomp_byte(compressed_values, (int)count, decompressed_values_byte, tilesize, blocksize);
    dbytes = (char *)decompressed_values_byte;
  } else if (bytepix == 2) {
    decompressed_values_short = (unsigned short *)malloc(tilesize * 2);
    fits_rdecomp_short(compressed_values, (int)count, decompressed_values_short, tilesize, blocksize);
    dbytes = (char *)decompressed_values_short;
  } else {
    decompressed_values_int = (unsigned int *)malloc(tilesize * 4);
    fits_rdecomp(compressed_values, (int)count, decompressed_values_int, tilesize, blocksize);
    dbytes = (char *)decompressed_values_int;
  }

  Py_END_ALLOW_THREADS

  if (PyErr_Occurred() != NULL) {
    // If an error condition inside the cfitsio function, the call inside
    // cfitsio should have called the ffpmsg function which sets the Python
    // exception, so we just return here to raise an error.
    return (PyObject *)NULL;
  }

  result = Py_BuildValue("y#", dbytes, tilesize * bytepix);
  free(dbytes);
  return result;
}

/* HCompress compression */

static PyObject *compress_hcompress_1_c(PyObject *self, PyObject *args) {

  const char *str;
  Py_ssize_t count;
  PyObject *result;

  int bytepix, nx, ny, scale;
  int status=0;  // Important to initialize this to zero otherwise will fail silently

  int maxelem;
  char *compressed_values;
  int *decompressed_values_int;
  long buffer_size;
  long long *decompressed_values_longlong;

  if (!PyArg_ParseTuple(args, "y#iiii", &str, &count, &nx, &ny, &scale, &bytepix)) {
    return NULL;
  }

  if (bytepix != 4 && bytepix != 8) {
    PyErr_SetString(PyExc_ValueError,
                    "HCompress can only work with 4 or 8 byte integers.");
    return (PyObject *)NULL;

  }

  if ((nx < 4) || (ny < 4)) {
    PyErr_SetString(PyExc_ValueError,
                    "HCOMPRESS requires tiles of at least 4x4 pixels.");
    return (PyObject *)NULL;
  }

  if (count != nx * ny * bytepix) {
    PyErr_SetString(PyExc_ValueError,
                    "The tile dimensions and dtype do not match the number of bytes provided.");
    return (PyObject *)NULL;
  }

  Py_BEGIN_ALLOW_THREADS

  // maxelem adapted from cfitsio's imcomp_calc_max_elem function
  maxelem = count / 4 * 2.2 + 26;

  // Apparently with the above calculation we can still end up allocating too
  // small of a buffer, this could never happen by more than 32 bytes
  // riiiiiight.
  // TODO: Do a small buffer calculation to tune this number like we did for PLIO
  compressed_values = (char *)calloc(maxelem + 4, sizeof(long long));
  buffer_size = (maxelem + 4) * sizeof(long long);

  if (bytepix == 4) {
    decompressed_values_int = (int *)str;
    fits_hcompress(decompressed_values_int, ny, nx, scale, compressed_values, &buffer_size, &status);
  } else {
    decompressed_values_longlong = (long long *)str;
    fits_hcompress64(decompressed_values_longlong, ny, nx, scale, compressed_values, &buffer_size, &status);
  }

  Py_END_ALLOW_THREADS

  if (PyErr_Occurred() != NULL) {
    // If an error condition inside the cfitsio function, the call inside
    // cfitsio should have called the ffpmsg function which sets the Python
    // exception, so we just return here to raise an error.
    return (PyObject *)NULL;
  }

  if (status != 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Status returned from cfitsio is not zero for an unknown reason.");
    return (PyObject *)NULL;
  }

  result = Py_BuildValue("y#", compressed_values, buffer_size);
  free(compressed_values);
  return result;
}

static PyObject *decompress_hcompress_1_c(PyObject *self, PyObject *args) {

  const unsigned char *str;
  char *dbytes;
  Py_ssize_t count;
  PyObject *result;

  int bytepix, nx, ny, scale, smooth;
  int status=0;  // Important to initialize this to zero otherwise will fail silently

  unsigned char *compressed_values;
  int *decompressed_values_int;
  long long *decompressed_values_longlong;

  if (!PyArg_ParseTuple(args, "y#iiiii", &str, &count, &nx, &ny, &scale, &smooth, &bytepix)) {
    return NULL;
  }

  if (bytepix != 4 && bytepix != 8) {
    PyErr_SetString(PyExc_ValueError,
                    "HCompress can only work with 4 or 8 byte integers.");
    return (PyObject *)NULL;

  }

  Py_BEGIN_ALLOW_THREADS

  compressed_values = (unsigned char *)str;

  dbytes = malloc(nx * ny * bytepix);

  if (bytepix == 4) {
    decompressed_values_int = (int *)dbytes;
    fits_hdecompress(compressed_values, smooth, decompressed_values_int, &ny, &nx, &scale, &status);
  } else {
    decompressed_values_longlong = (long long *)dbytes;
    fits_hdecompress64(compressed_values, smooth, decompressed_values_longlong, &ny, &nx, &scale, &status);
  }

  Py_END_ALLOW_THREADS

  if (PyErr_Occurred() != NULL) {
    // If an error condition inside the cfitsio function, the call inside
    // cfitsio should have called the ffpmsg function which sets the Python
    // exception, so we just return here to raise an error.
    return (PyObject *)NULL;
  }

  if (status != 0) {
    PyErr_SetString(PyExc_ValueError,
                    "Status returned from cfitsio is not zero for an unknown reason.");
    return (PyObject *)NULL;
  }

  // fits_hdecompress[64] always returns 4 byte integers
  result = Py_BuildValue("y#", dbytes, nx * ny * 4);
  free(dbytes);
  return result;
}

static PyObject *quantize_float_c(PyObject *self, PyObject *args) {

  const char *input_bytes;
  Py_ssize_t nbytes;
  PyObject *result;

  float *input_data;

  long row, nx, ny;
  int nullcheck;
  float in_null_value;
  float qlevel;
  int dither_method;

  int *quantized_data;
  char *quantized_bytes;
  double bscale, bzero;
  int iminval, imaxval;

  int status;

  Py_ssize_t output_length;

  if (!PyArg_ParseTuple(args, "y#lllidfi", &input_bytes, &nbytes, &row, &nx,
                        &ny, &nullcheck, &in_null_value, &qlevel,
                        &dither_method)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS

  input_data = (float *)input_bytes;
  quantized_data = (int *)malloc(nx * ny * sizeof(int));

  status = fits_quantize_float(row, input_data, nx, ny, nullcheck, in_null_value, qlevel,
                               dither_method, quantized_data, &bscale, &bzero, &iminval,
                               &imaxval);

  quantized_bytes = (char *)quantized_data;

  output_length = nx * ny * sizeof(int);

  Py_END_ALLOW_THREADS

  result = Py_BuildValue("y#iddii", quantized_bytes, output_length, status,
                         bscale, bzero, iminval, imaxval);
  free(quantized_bytes);
  return result;
}

static PyObject *quantize_double_c(PyObject *self, PyObject *args) {

  const char *input_bytes;
  Py_ssize_t nbytes;
  PyObject *result;

  double *input_data;

  long row, nx, ny;
  int nullcheck;
  double in_null_value;
  float qlevel;
  int dither_method;

  int *quantized_data;
  char *quantized_bytes;
  double bscale, bzero;
  int iminval, imaxval;

  int status;

  if (!PyArg_ParseTuple(args, "y#lllidfi", &input_bytes, &nbytes, &row, &nx,
                        &ny, &nullcheck, &in_null_value, &qlevel,
                        &dither_method)) {
    return NULL;
  }

  Py_BEGIN_ALLOW_THREADS

  input_data = (double *)input_bytes;
  quantized_data = (int *)malloc(nx * ny * sizeof(int));

  status = fits_quantize_double(row, input_data, nx, ny, nullcheck, in_null_value,
                                qlevel, dither_method, quantized_data, &bscale, &bzero,
                                &iminval, &imaxval);

  quantized_bytes = (char *)quantized_data;

  Py_END_ALLOW_THREADS

  result = Py_BuildValue("y#iddii", quantized_bytes, nx * ny * sizeof(int), status,
                                    bscale, bzero, iminval, imaxval);
  free(quantized_bytes);
  return result;
}

static PyObject *unquantize_float_c(PyObject *self, PyObject *args) {

  const char *input_bytes;
  Py_ssize_t nbytes;
  PyObject *result;

  long row, npix;
  int nullcheck;
  int tnull;
  float nullval;
  int dither_method;

  double bscale, bzero;
  int bytepix; // int size
  int status = 0;

  int *anynull;
  float *output_data;
  char *output_bytes;

  if (!PyArg_ParseTuple(args, "y#llddiiifi", &input_bytes, &nbytes, &row, &npix,
                        &bscale, &bzero, &dither_method, &nullcheck, &tnull,
                        &nullval, &bytepix)) {
    return NULL;
  }

  // TODO: add support, if needed, for nullcheck=1

  Py_BEGIN_ALLOW_THREADS

  anynull = (int *)malloc(npix * sizeof(int));
  output_data = (float *)calloc(npix, sizeof(float));

  if (bytepix == 1) {
      unquantize_i1r4(row, (unsigned char *)input_bytes, npix, bscale, bzero,
                      dither_method, nullcheck, (unsigned char)tnull, nullval,
                      NULL, anynull, output_data, &status);
  } else if (bytepix == 2) {
      unquantize_i2r4(row, (short *)input_bytes, npix, bscale, bzero,
                      dither_method, nullcheck, (short)tnull, nullval, NULL,
                      anynull, output_data, &status);
  } else if (bytepix == 4) {
      unquantize_i4r4(row, (int *)input_bytes, npix, bscale, bzero, dither_method,
                      nullcheck, (int)tnull, nullval, NULL, anynull, output_data,
                      &status);
  }

  output_bytes = (char *)output_data;

  Py_END_ALLOW_THREADS

  result = Py_BuildValue("y#", output_bytes, npix * sizeof(float));
  free(output_bytes);
  free(anynull);
  return result;
}

static PyObject *unquantize_double_c(PyObject *self, PyObject *args) {

  const char *input_bytes;
  Py_ssize_t nbytes;
  PyObject *result;

  long row, npix;
  int nullcheck;
  int tnull;
  double nullval;
  int dither_method;

  double bscale, bzero;
  int bytepix; // int size
  int status = 0;

  int *anynull;
  double *output_data;
  char *output_bytes;

  if (!PyArg_ParseTuple(args, "y#llddiiidi", &input_bytes, &nbytes, &row, &npix,
                        &bscale, &bzero, &dither_method, &nullcheck, &tnull,
                        &nullval, &bytepix)) {
    return NULL;
  }

  // TODO: add support, if needed, for nullcheck=1

  Py_BEGIN_ALLOW_THREADS

  anynull = (int *)malloc(npix * sizeof(int));
  output_data = (double *)malloc(npix * sizeof(double));

  if (bytepix == 1) {
      unquantize_i1r8(row, (unsigned char *)input_bytes, npix, bscale, bzero,
                      dither_method, nullcheck, (unsigned char)tnull, nullval,
                      NULL, anynull, output_data, &status);
  } else if (bytepix == 2) {
      unquantize_i2r8(row, (short *)input_bytes, npix, bscale, bzero,
                      dither_method, nullcheck, (short)tnull, nullval, NULL,
                      anynull, output_data, &status);
  } else if (bytepix == 4) {
      unquantize_i4r8(row, (int *)input_bytes, npix, bscale, bzero, dither_method,
                      nullcheck, (int)tnull, nullval, NULL, anynull, output_data,
                      &status);
  }

  output_bytes = (char *)output_data;

  Py_END_ALLOW_THREADS

  result = Py_BuildValue("y#", output_bytes, npix * sizeof(double));
  free(output_bytes);
  free(anynull);
  return result;
}
