/******************************************************************************
 * C extension code for astropy.utils.xml.iterparse
 *
 * Everything in this file has an alternate Python implementation and
 * is included for performance reasons only.
 *
 * It has two main parts:
 *
 *   - An IterParser object which parses an XML file using the expat
 *     library, feeding expat events through a Python iterator.  It is
 *     faster and more memory efficient than the alternatives in the
 *     Python standard library because it does not build a tree of
 *     objects, and also throws away most text nodes, since for
 *     astropy.io.votable (the primary user of this library) we only
 *     care about simple text nodes contained between a single pair of
 *     open/close element nodes.  It also has an optimization for
 *     recognizing the most commonly occurring element in a VO file,
 *     "TD".
 *
 *   - Two functions, escape_xml() and escape_xml_cdata() that escape
 *     XML much faster than the alternatives in the Python standard
 *     library.
 ******************************************************************************/

#include <stdio.h>
#include <Python.h>
#include "structmember.h"

#include "expat.h"

/******************************************************************************
 * Convenience macros and functions
 ******************************************************************************/
#ifdef _MSC_VER
#define inline
#endif

#undef  CLAMP
#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

static Py_ssize_t
next_power_of_2(Py_ssize_t n)
{
    /* Calculate the next-higher power of two that is >= 'n' */

    /* These instructions are intended for uint32_t and originally
       from http://www-graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
       Py_ssize_t is the same as the C datatype ssize_t and is
       32/64-bits on 32-bit and 64-bit systems respectively. The implementation
       here accounts for both.

       Limitations: Since a signed size_t (ssize_t) was required for the underlying CPython implementation,
       on 32 bit systems, it will not be possible to allocate memory sizes between
       SSIZE_MAX and SIZE_MAX (i.e, for allocations in the range [2^31, 2^32 - 1] bytes)
       even though such memory might be available and accessible on the (32-bit) computer. That
       said, since the underlying CPython implementation *also* uses Py_ssize_t (i.e., ssize_t),
       it is safe to assume that such memory allocations would probably not be usable anyway.

       TLDR: Will work on 64-bit machines but be careful on 32 bit machines when reading in
       ~ 2+ GB of memory -- @manodeep 2020-03-27
    */

    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    if(sizeof(Py_ssize_t) > 4) {
        n |= n >> 32; /* this works for 64-bit systems but will need to be updated if (Py)_ssize_t
                         ever increases beyond 64-bits */
    }
    n++;

    return n;
}

/******************************************************************************
 * Python version compatibility macros
 ******************************************************************************/

#if BYTEORDER == 1234
# define TD_AS_INT      0x00004454
# define TD_AS_INT_MASK 0x00ffffff
#else
# define TD_AS_INT      0x54440000
# define TD_AS_INT_MASK 0xffffff00
#endif

/* Clang doesn't like the hackish stuff PyTuple_SET_ITEM does... */
#ifdef __clang__
#undef PyTuple_SET_ITEM
#define PyTuple_SET_ITEM(a, b, c) PyTuple_SetItem((a), (b), (c))
#endif

/******************************************************************************
 * IterParser type
 ******************************************************************************/
typedef struct {
    PyObject_HEAD
    XML_Parser parser;          /* The expat parser */
    int        done;            /* True when expat parser has read to EOF */

    /* File-like object reading */
    PyObject*  fd;              /* Python file object */
    int        file;            /* C file descriptor */
    PyObject*  read;            /* The read method on the file object */
    Py_ssize_t    buffersize;      /* The size of the read buffer */
    XML_Char*  buffer;          /* The read buffer */

    /* Text nodes */
    Py_ssize_t text_alloc;      /* The allocated size of the text buffer */
    Py_ssize_t text_size;       /* The size of the content in the text buffer */
    XML_Char*  text;            /* Text buffer (for returning text nodes) */
    int        keep_text;       /* Flag: keep appending text chunks to the current text node */

    /* XML event queue */
    PyObject** queue;
    Py_ssize_t queue_size;
    Py_ssize_t queue_read_idx;
    Py_ssize_t queue_write_idx;

    /* Store the last Python exception so it can be returned when
       dequeuing events */
    PyObject*  error_type;
    PyObject*  error_value;
    PyObject*  error_traceback;

    /* Store the position for any XML exceptions that may be
       returned later */
    unsigned long last_line;
    unsigned long last_col;

    /* "Constants" for efficiency */
    PyObject*  dict_singleton;  /* Empty dict */
    PyObject*  td_singleton;    /* String "TD" */
    PyObject*  read_args;       /* (buffersize) */
} IterParser;

/******************************************************************************
 * Tuple queue
 ******************************************************************************/

/**
 * Extend the tuple queue based on the new length of the textual XML input.
 * This helps to cope with situations where the input is longer than
 * requested (as occurs with transparent decompression of the input
 * stream), and for the initial allocation to combine the logic in one place.
 */
static int
queue_realloc(IterParser *self, Py_ssize_t req_size)
{
    PyObject** new_queue;
    Py_ssize_t n = req_size / 2;

    if (n <= self->queue_size)
        return 0;

    new_queue = realloc(self->queue, sizeof(PyObject*) * (size_t)n);

    if (new_queue == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory for XML parsing queue.");
        /*
         * queue_realloc() is only called from IterParser_init() or
         * IterParser_next() in situations where the queue is clear
         * and empty.  If this function were to be used in other
         * situations it would be wise to iterate over the queue and
         * clear/decrement the individual references, to save work for
         * the garbage collector (in an out-of-memory situation).
         */
        goto fail;
    }

    self->queue = new_queue;
    self->queue_size = n;
    return 0;

fail:
    free(self->queue);
    self->queue = NULL;
    self->queue_size = 0;
    return -1;
}

/******************************************************************************
 * Text buffer
 ******************************************************************************/

/**
 * Reallocate text buffer to the next highest power of two that fits the
 * requested size.
 */
static int
text_realloc(IterParser *self, Py_ssize_t req_size)
{
    Py_ssize_t  n       = req_size;
    char       *new_mem = NULL;

    if (req_size < self->text_alloc) {
        return 0;
    }

    /* Calculate the next-highest power of two */
    n = next_power_of_2(n);

    if (n < req_size) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory for XML text.");
        return -1;
    }

    new_mem = malloc(n * sizeof(XML_Char));
    if (new_mem == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory for XML text.");
        return -1;
    }

    memcpy(new_mem, self->text, (size_t)(self->text_size + 1) * sizeof(XML_Char));

    free(self->text);
    self->text = new_mem;
    self->text_alloc = n;

    return 0;
}

#define IS_WHITESPACE(c) ((c) == (XML_Char)0x20 || \
                          (c) == (XML_Char)0x0d || \
                          (c) == (XML_Char)0x0a || \
                          (c) == (XML_Char)0x09)

/*
 * Append text to the text buffer.
 *
 * For the first chunk of text, all whitespace characters before the
 * first non-whitespace character are stripped.  This saves time
 * stripping on the Python side later.
 */
static int
text_append(IterParser *self, const XML_Char *data, Py_ssize_t len)
{
    Py_ssize_t new_size;

    if (len == 0) {
        return 0;
    }

    /* If this is the first chunk, handle whitespace */
    if (self->text_size == 0) {
        while (len && IS_WHITESPACE(*data)) {
            ++data;
            --len;
        }
    }

    /* Grow text buffer if necessary */
    new_size = self->text_size + len;
    if (text_realloc(self, new_size + 1)) {
        return -1;
    }

    memcpy(self->text + self->text_size,
           data,
           (size_t)len * sizeof(XML_Char));

    self->text_size = new_size;
    self->text[self->text_size] = (XML_Char)0x0;

    return 0;
}

/*
 * Erase all content from the text buffer.
 */
static void
text_clear(IterParser *self)
{
    self->text[0] = (XML_Char)0;
    self->text_size = 0;
}

/******************************************************************************
 * XML event handling
 ******************************************************************************/

/*
 * Make a "position tuple" from the current expat parser state.  This
 * is used to communicate the position of the parser within the file
 * to the Python side for generation of meaningful error messages.
 *
 * It is of the form (line, col), where line and col are both PyInts.
 */
static inline PyObject*
make_pos(const IterParser *self)
{
    return Py_BuildValue(
            "(nn)",
            (size_t)self->last_line,
            (size_t)self->last_col);
}

/*
 * Removes the namespace from an element or attribute name, that is,
 * remove everything before the first colon.  The namespace is not
 * needed to parse standards-compliant VOTable files.
 *
 * The returned pointer is an internal pointer to the buffer passed
 * in.
 */
static const XML_Char *
remove_namespace(const XML_Char *name)
{
    const XML_Char*  name_start = NULL;

    /* If there is a namespace specifier, just chop it off */
    for (name_start = name; *name_start != '\0'; ++name_start) {
        if (*name_start == ':') {
            break;
        }
    }

    if (*name_start == ':') {
        ++name_start;
    } else {
        name_start = name;
    }

    return name_start;
}

/*
 * Handle the expat startElement event.
 */
static void
startElement(IterParser *self, const XML_Char *name, const XML_Char **atts)
{
    PyObject*        pyname = NULL;
    PyObject*        pyatts = NULL;
    const XML_Char** att_ptr = atts;
    const XML_Char*  name_start = NULL;
    PyObject*        tuple = NULL;
    PyObject*        key = NULL;
    PyObject*        val = NULL;
    PyObject*        pos = NULL;

    /* If we've already had an error in a previous call, don't make
       things worse. */
    if (PyErr_Occurred() != NULL) {
        XML_StopParser(self->parser, 0);
        return;
    }

    /* Don't overflow the queue -- in practice this should *never* happen */
    if (self->queue_write_idx < self->queue_size) {
        tuple = PyTuple_New(4);
        if (tuple == NULL) {
            goto fail;
        }

        Py_INCREF(Py_True);
        PyTuple_SET_ITEM(tuple, 0, Py_True);

        /* This is an egregious but effective optimization.  Since by
           far the most frequently occurring element name in a large
           VOTABLE file is TD, we explicitly check for it here with
           integer comparison to avoid the lookup in the interned
           string table in PyString_InternFromString, and return a
           singleton string for "TD" */
        if ((*(int*)name & TD_AS_INT_MASK) == TD_AS_INT) {
            Py_INCREF(self->td_singleton);
            PyTuple_SetItem(tuple, 1, self->td_singleton);
        } else {
            name_start = remove_namespace(name);

            pyname = PyUnicode_FromString(name_start);
            if (pyname == NULL) {
                goto fail;
            }
            PyTuple_SetItem(tuple, 1, pyname);
            pyname = NULL;
        }

        if (*att_ptr) {
            pyatts = PyDict_New();
            if (pyatts == NULL) {
                goto fail;
            }
            do {
                key = PyUnicode_FromString(*att_ptr);
                if (key == NULL) {
                    goto fail;
                }
                val = PyUnicode_FromString(*(att_ptr + 1));
                if (val == NULL) {
                    Py_DECREF(key);
                    goto fail;
                }
                if (PyDict_SetItem(pyatts, key, val)) {
                    Py_DECREF(key);
                    Py_DECREF(val);
                    goto fail;
                }
                Py_DECREF(key);
                Py_DECREF(val);
                key = val = NULL;

                att_ptr += 2;
            } while (*att_ptr);
        } else {
            Py_INCREF(self->dict_singleton);
            pyatts = self->dict_singleton;
        }

        PyTuple_SetItem(tuple, 2, pyatts);
        pyatts = NULL;

        self->last_line = (unsigned long)XML_GetCurrentLineNumber(
            self->parser);
        self->last_col = (unsigned long)XML_GetCurrentColumnNumber(
            self->parser);

        pos = make_pos(self);
        if (pos == NULL) {
            goto fail;
        }
        PyTuple_SetItem(tuple, 3, pos);
        pos = NULL;

        text_clear(self);

        self->keep_text = 1;

        self->queue[self->queue_write_idx++] = tuple;
    } else {
        PyErr_SetString(
            PyExc_RuntimeError,
            "XML queue overflow in startElement.  This most likely indicates an internal bug.");
        goto fail;
    }

    return;

 fail:
    Py_XDECREF(tuple);
    Py_XDECREF(pyatts);
    XML_StopParser(self->parser, 0);
}

/*
 * Handle the expat endElement event.
 */
static void
endElement(IterParser *self, const XML_Char *name)
{
    PyObject*       pyname     = NULL;
    PyObject*       tuple      = NULL;
    PyObject*       pytext     = NULL;
    const XML_Char* name_start = NULL;
    XML_Char*       end;
    PyObject*       pos        = NULL;

    /* If we've already had an error in a previous call, don't make
       things worse. */
    if (PyErr_Occurred() != NULL) {
        XML_StopParser(self->parser, 0);
        return;
    }

    /* Don't overflow the queue -- in practice this should *never* happen */
    if (self->queue_write_idx < self->queue_size) {
        tuple = PyTuple_New(4);
        if (tuple == NULL) {
            goto fail;
        }

        Py_INCREF(Py_False);
        PyTuple_SET_ITEM(tuple, 0, Py_False);

        /* This is an egregious but effective optimization.  Since by
           far the most frequently occurring element name in a large
           VOTABLE file is TD, we explicitly check for it here with
           integer comparison to avoid the lookup in the interned
           string table in PyString_InternFromString, and return a
           singleton string for "TD" */
        if ((*(int*)name & TD_AS_INT_MASK) == TD_AS_INT) {
            Py_INCREF(self->td_singleton);
            PyTuple_SetItem(tuple, 1, self->td_singleton);
        } else {
            name_start = remove_namespace(name);

            pyname = PyUnicode_FromString(name_start);
            if (pyname == NULL) {
                goto fail;
            }
            PyTuple_SetItem(tuple, 1, pyname);
            pyname = NULL;
        }

        /* Cut whitespace off the end of the string */
        end = self->text + self->text_size - 1;
        while (end >= self->text && IS_WHITESPACE(*end)) {
            --end;
            --self->text_size;
        }

        pytext = PyUnicode_FromStringAndSize(self->text, self->text_size);
        if (pytext == NULL) {
            goto fail;
        }
        PyTuple_SetItem(tuple, 2, pytext);
        pytext = NULL;

        pos = make_pos(self);
        if (pos == NULL) {
            goto fail;
        }
        PyTuple_SetItem(tuple, 3, pos);
        pos = NULL;

        self->keep_text = 0;

        self->queue[self->queue_write_idx++] = tuple;
    } else {
        PyErr_SetString(
            PyExc_RuntimeError,
            "XML queue overflow in endElement.  This most likely indicates an internal bug.");
        goto fail;
    }

    return;

 fail:
    Py_XDECREF(tuple);
    XML_StopParser(self->parser, 0);
}

/*
 * Handle the expat characterData event.
 */
static void
characterData(IterParser *self, const XML_Char *text, int len)
{
    /* If we've already had an error in a previous call, don't make
       things worse. */
    if (PyErr_Occurred() != NULL) {
        XML_StopParser(self->parser, 0);
        return;
    }

    if (self->text_size == 0) {
        self->last_line = (unsigned long)XML_GetCurrentLineNumber(
            self->parser);
        self->last_col = (unsigned long)XML_GetCurrentColumnNumber(
            self->parser);
    }

    if (self->keep_text) {
        (void)text_append(self, text, (Py_ssize_t)len);
    }
}

/*
 * Handle the XML declaration so that we can determine its encoding.
 */
static void
xmlDecl(IterParser *self, const XML_Char *version,
        const XML_Char *encoding, int standalone)
{
    PyObject* tuple        = NULL;
    PyObject* xml_str      = NULL;
    PyObject* attrs        = NULL;
    PyObject* encoding_str = NULL;
    PyObject* version_str  = NULL;
    PyObject* pos          = NULL;

    if (self->queue_write_idx < self->queue_size) {
        tuple = PyTuple_New(4);
        if (tuple == NULL) {
            goto fail;
        }

        Py_INCREF(Py_True);
        PyTuple_SET_ITEM(tuple, 0, Py_True);

        xml_str = PyUnicode_FromString("xml");
        if (xml_str == NULL) {
            goto fail;
        }
        PyTuple_SET_ITEM(tuple, 1, xml_str);
        xml_str = NULL;

        attrs = PyDict_New();
        if (attrs == NULL) {
            goto fail;
        }

        if (encoding) {
            encoding_str = PyUnicode_FromString(encoding);
        } else {
            encoding_str = PyUnicode_FromString("");
        }
        if (encoding_str == NULL) {
            goto fail;
        }
        if (PyDict_SetItemString(attrs, "encoding", encoding_str)) {
            Py_DECREF(encoding_str);
            goto fail;
        }
        Py_DECREF(encoding_str);
        encoding_str = NULL;

        if (version) {
            version_str = PyUnicode_FromString(version);
        } else {
            version_str = PyUnicode_FromString("");
        }
        if (version_str == NULL) {
            goto fail;
        }
        if (PyDict_SetItemString(attrs, "version", version_str)) {
            Py_DECREF(version_str);
            goto fail;
        }
        Py_DECREF(version_str);
        version_str = NULL;

        PyTuple_SET_ITEM(tuple, 2, attrs);
        attrs = NULL;

        self->last_line = (unsigned long)XML_GetCurrentLineNumber(
            self->parser);
        self->last_col = (unsigned long)XML_GetCurrentColumnNumber(
            self->parser);

        pos = make_pos(self);
        if (pos == NULL) {
            goto fail;
        }
        PyTuple_SetItem(tuple, 3, pos);
        pos = NULL;

        self->queue[self->queue_write_idx++] = tuple;
    } else {
        PyErr_SetString(
            PyExc_RuntimeError,
            "XML queue overflow in xmlDecl.  This most likely indicates an internal bug.");
        goto fail;
    }

    return;

 fail:
    Py_XDECREF(tuple);
    Py_XDECREF(attrs);
    XML_StopParser(self->parser, 0);
}

/*
 * The object itself is an iterator, just return self for "iter(self)"
 * on the Python side.
 */
static PyObject *
IterParser_iter(IterParser* self)
{
    Py_INCREF(self);
    return (PyObject*) self;
}

/*
 * Get the next element from the iterator.
 *
 * The expat event handlers above (startElement, endElement, characterData) add
 * elements to the queue, which are then dequeued by this method.
 *
 * Care must be taken to store and later raise exceptions.  Any
 * exceptions raised in the expat callbacks must be stored and then
 * later thrown once the queue is emptied, otherwise the exception is
 * raised "too early" in queue order.
 */
static PyObject *
IterParser_next(IterParser* self)
{
    PyObject*  data = NULL;
    XML_Char*  buf;
    Py_ssize_t buflen;

    /* Is there anything in the queue to return? */
    if (self->queue_read_idx < self->queue_write_idx) {
        return self->queue[self->queue_read_idx++];
    }

    /* Now that the queue is empty, is there an error we need to raise? */
    if (self->error_type) {
        PyErr_Restore(self->error_type, self->error_value, self->error_traceback);
        self->error_type = NULL;
        self->error_value = NULL;
        self->error_traceback = NULL;
        return NULL;
    }

    /* The queue is empty -- have we already fed the entire file to
       expat?  If so, we are done and indicate the end of the iterator
       by simply returning NULL. */
    if (self->done) {
        return NULL;
    }

    self->queue_read_idx = 0;
    self->queue_write_idx = 0;

    do {
        /* Handle a generic Python read method */
        if (self->read) {
            data = PyObject_CallObject(self->read, self->read_args);
            if (data == NULL) {
                goto fail;
            }

            if (PyBytes_AsStringAndSize(data, &buf, &buflen) == -1) {
                Py_DECREF(data);
                goto fail;
            }

            if (buflen < self->buffersize) {
                /* EOF detection method only works for local regular files */
                self->done = 1;
            }
        /* Handle a real C file descriptor or handle -- this is faster
           if we've got one. */
        } else {
            buflen = (Py_ssize_t)read(
                self->file, self->buffer, (size_t)self->buffersize);
            if (buflen == -1) {
                PyErr_SetFromErrno(PyExc_OSError);
                goto fail;
            } else if (buflen < self->buffersize) {
                /* EOF detection method only works for local regular files */
                self->done = 1;
            }

            buf = self->buffer;
        }

        if(queue_realloc(self, buflen)) {
            Py_XDECREF(data);
            goto fail;
        }

        /* Feed the read buffer to expat, which will call the event handlers */
        if (XML_Parse(self->parser, buf, (int)buflen, self->done) == XML_STATUS_ERROR) {
            /* One of the event handlers raised a Python error, make
               note of it -- it won't be thrown until the queue is
               emptied. */
            if (PyErr_Occurred() != NULL) {
                goto fail;
            }

            /* expat raised an error, make note of it -- it won't be thrown
               until the queue is emptied. */
            Py_XDECREF(data);
            PyErr_Format(
                PyExc_ValueError, "%lu:%lu: %s",
                XML_GetCurrentLineNumber(self->parser),
                XML_GetCurrentColumnNumber(self->parser),
                XML_ErrorString(XML_GetErrorCode(self->parser)));
            goto fail;
        }
        Py_XDECREF(data);

        if (PyErr_Occurred() != NULL) {
            goto fail;
        }
    } while (self->queue_write_idx == 0 && self->done == 0);

    if (self->queue_write_idx == 0) {
        return NULL;
    }

    if (self->queue_write_idx >= self->queue_size) {
        PyErr_SetString(
            PyExc_RuntimeError,
            "XML queue overflow.  This most likely indicates an internal bug.");
        return NULL;
    }

    return self->queue[self->queue_read_idx++];

 fail:
    /* We got an exception somewhere along the way.  Store the exception in
       the IterParser object, but clear the exception in the Python interpreter,
       so we can empty the event queue and raise the exception later. */
    PyErr_Fetch(&self->error_type, &self->error_value, &self->error_traceback);
    PyErr_Clear();

    if (self->queue_read_idx < self->queue_write_idx) {
        return self->queue[self->queue_read_idx++];
    }

    PyErr_Restore(self->error_type, self->error_value, self->error_traceback);
    self->error_type = NULL;
    self->error_value = NULL;
    self->error_traceback = NULL;
    return NULL;
}

/******************************************************************************
 * IterParser object lifetime
 ******************************************************************************/

/* To support cyclical garbage collection, all PyObject's must be
   visited. */
static int
IterParser_traverse(IterParser *self, visitproc visit, void *arg)
{
    int vret;
    Py_ssize_t read_index;

    read_index = self->queue_read_idx;
    while (read_index < self->queue_write_idx) {
        vret = visit(self->queue[read_index++], arg);
        if (vret != 0) return vret;
    }

    if (self->fd) {
        vret = visit(self->fd, arg);
        if (vret != 0) return vret;
    }

    if (self->read) {
        vret = visit(self->read, arg);
        if (vret != 0) return vret;
    }

    if (self->read_args) {
        vret = visit(self->read_args, arg);
        if (vret != 0) return vret;
    }

    if (self->dict_singleton) {
        vret = visit(self->dict_singleton, arg);
        if (vret != 0) return vret;
    }

    if (self->td_singleton) {
        vret = visit(self->td_singleton, arg);
        if (vret != 0) return vret;
    }

    if (self->error_type) {
        vret = visit(self->error_type, arg);
        if (vret != 0) return vret;
    }

    if (self->error_value) {
        vret = visit(self->error_value, arg);
        if (vret != 0) return vret;
    }

    if (self->error_traceback) {
        vret = visit(self->error_traceback, arg);
        if (vret != 0) return vret;
    }

    return 0;
}

/* To support cyclical garbage collection */
static int
IterParser_clear(IterParser *self)
{
    PyObject *tmp;

    while (self->queue_read_idx < self->queue_write_idx) {
        tmp = self->queue[self->queue_read_idx];
        self->queue[self->queue_read_idx] = NULL;
        Py_XDECREF(tmp);
        self->queue_read_idx++;
    }

    tmp = self->fd;
    self->fd = NULL;
    Py_XDECREF(tmp);

    tmp = self->read;
    self->read = NULL;
    Py_XDECREF(tmp);

    tmp = self->read_args;
    self->read_args = NULL;
    Py_XDECREF(tmp);

    tmp = self->dict_singleton;
    self->dict_singleton = NULL;
    Py_XDECREF(tmp);

    tmp = self->td_singleton;
    self->td_singleton = NULL;
    Py_XDECREF(tmp);

    tmp = self->error_type;
    self->error_type = NULL;
    Py_XDECREF(tmp);

    tmp = self->error_value;
    self->error_value = NULL;
    Py_XDECREF(tmp);

    tmp = self->error_traceback;
    self->error_traceback = NULL;
    Py_XDECREF(tmp);

    return 0;
}

/*
 * Deallocate the IterParser object.  For the internal PyObject*, just
 * punt to IterParser_clear.
 */
static void
IterParser_dealloc(IterParser* self)
{
    IterParser_clear(self);

    free(self->buffer); self->buffer = NULL;
    free(self->queue);  self->queue = NULL;
    free(self->text);   self->text = NULL;
    if (self->parser != NULL) {
        XML_ParserFree(self->parser);
        self->parser = NULL;
    }

    Py_TYPE(self)->tp_free((PyObject*)self);
}

/*
 * Initialize the memory for an IterParser object
 */

static PyObject *
IterParser_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    IterParser *self = NULL;

    self = (IterParser *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->parser          = NULL;
        self->fd              = NULL;
        self->file            = -1;
        self->read            = NULL;
        self->read_args       = NULL;
        self->dict_singleton  = NULL;
        self->td_singleton    = NULL;
        self->buffersize      = 0;
        self->buffer          = NULL;
        self->queue_read_idx  = 0;
        self->queue_write_idx = 0;
        self->text_alloc      = 0;
        self->text_size       = 0;
        self->text            = NULL;
        self->keep_text       = 0;
        self->done            = 0;
        self->queue_size      = 0;
        self->queue           = NULL;
        self->error_type      = NULL;
        self->error_value     = NULL;
        self->error_traceback = NULL;
    }

    return (PyObject *)self;
}

/*
 * Initialize an IterParser object
 *
 * The Python arguments are:
 *
 *    *fd*: A Python file object or a callable object
 *    *buffersize*: The size of the read buffer
 */
static int
IterParser_init(IterParser *self, PyObject *args, PyObject *kwds)
{
    PyObject* fd              = NULL;
    PyObject* read            = NULL;
    Py_ssize_t   buffersize      = 1 << 14;

    static char *kwlist[] = {"fd", "buffersize", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|n:IterParser.__init__", kwlist,
                                     &fd, &buffersize)) {
        return -1;
    }

    /* Keep the buffersize within a reasonable range */
    self->buffersize = CLAMP(buffersize, (Py_ssize_t)(1 << 10), (Py_ssize_t)(1 << 24));
#ifdef __clang__
    /* Clang can't handle the file descriptors Python gives us,
       so in that case, we just call the object's read method. */
    read = PyObject_GetAttrString(fd, "read");
    if (read != NULL) {
        fd = read;
    }
#else
    self->file = PyObject_AsFileDescriptor(fd);
    if (self->file != -1) {
        /* This is a real C file handle or descriptor.  We therefore
           need to allocate our own read buffer, and get the real C
           object. */
        self->buffer = malloc((size_t)self->buffersize);
        if (self->buffer == NULL) {
            PyErr_SetString(PyExc_MemoryError, "Out of memory");
            goto fail;
        }
        self->fd = fd;   Py_INCREF(self->fd);
        lseek(self->file, 0, SEEK_SET);
    } else
#endif
    if (PyCallable_Check(fd)) {
        /* fd is a Python callable */
        self->fd = fd;   Py_INCREF(self->fd);
        self->read = fd; Py_INCREF(self->read);
    } else {
        PyErr_SetString(
            PyExc_TypeError,
            "Arg 1 to iterparser must be a file object or callable object");
        goto fail;
    }

    PyErr_Clear();

    self->queue_read_idx  = 0;
    self->queue_write_idx = 0;
    self->done            = 0;

    self->text = malloc((size_t)buffersize * sizeof(XML_Char));
    self->text_alloc = buffersize;
    if (self->text == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        goto fail;
    }
    text_clear(self);

    self->read_args = Py_BuildValue("(n)", buffersize);
    if (self->read_args == NULL) {
        goto fail;
    }

    self->dict_singleton = PyDict_New();
    if (self->dict_singleton == NULL) {
        goto fail;
    }

    self->td_singleton = PyUnicode_FromString("TD");
    if (self->td_singleton == NULL) {
        goto fail;
    }

    if (queue_realloc(self, buffersize)) {
        goto fail;
    }

    /* Set up an expat parser with our callbacks */
    self->parser = XML_ParserCreate(NULL);
    if (self->parser == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        goto fail;
    }
    XML_SetUserData(self->parser, self);
    XML_SetElementHandler(
        self->parser,
        (XML_StartElementHandler)startElement,
        (XML_EndElementHandler)endElement);
    XML_SetCharacterDataHandler(
        self->parser,
        (XML_CharacterDataHandler)characterData);
    XML_SetXmlDeclHandler(
        self->parser,
        (XML_XmlDeclHandler)xmlDecl);

    Py_XDECREF(read);

    return 0;

 fail:
    Py_XDECREF(read);
    Py_XDECREF(self->fd);
    Py_XDECREF(self->read);
    free(self->text);
    Py_XDECREF(self->dict_singleton);
    Py_XDECREF(self->td_singleton);
    Py_XDECREF(self->read_args);
    free(self->queue);

    return -1;
}

static PyMemberDef IterParser_members[] =
{
    {NULL}  /* Sentinel */
};

static PyMethodDef IterParser_methods[] =
{
    {NULL}  /* Sentinel */
};

static PyTypeObject IterParserType =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    "astropy.utils.xml._iterparser.IterParser",    /*tp_name*/
    sizeof(IterParser),         /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)IterParser_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /*tp_flags*/
    "IterParser objects",       /* tp_doc */
    (traverseproc)IterParser_traverse, /* tp_traverse */
    (inquiry)IterParser_clear,  /* tp_clear */
    0,                          /* tp_richcompare */
    0,                          /* tp_weaklistoffset */
    (getiterfunc)IterParser_iter, /* tp_iter */
    (iternextfunc)IterParser_next, /* tp_iternext */
    IterParser_methods,         /* tp_methods */
    IterParser_members,         /* tp_members */
    0,                          /* tp_getset */
    0,                          /* tp_base */
    0,                          /* tp_dict */
    0,                          /* tp_descr_get */
    0,                          /* tp_descr_set */
    0,                          /* tp_dictoffset */
    (initproc)IterParser_init,  /* tp_init */
    0,                          /* tp_alloc */
    IterParser_new,             /* tp_new */
};

/******************************************************************************
 * XML escaping
 ******************************************************************************/

/* These are in reverse order by input character */
static const char* escapes_cdata[] = {
    ">", "&gt;",
    "<", "&lt;",
    "&", "&amp;",
    "\0", "\0",
};

/* These are in reverse order by input character */
static const char* escapes[] = {
    ">", "&gt;",
    "<", "&lt;",
    "'", "&apos;",
    "&", "&amp;",
    "\"", "&quot;",
    "\0", "\0"
};

/* Implementation of escape_xml.
 *
 * Returns:
 *  * 0  : No need to escape
 *  * >0 : output is escaped
 *  * -1 : error
 */
static Py_ssize_t
_escape_xml_impl(const char *input, Py_ssize_t input_len,
                 char **output, const char **escapes)
{
    Py_ssize_t i;
    int count = 0;
    char *p = NULL;
    const char** esc;
    const char* ent;

    for (i = 0; i < input_len; ++i) {
        for (esc = escapes; ; esc += 2) {
            if ((unsigned char)input[i] > **esc) {
                break;
            } else if (input[i] == **esc) {
                ++count;
                break;
            }
        }
    }

    if (!count) {
        return 0;
    }

    p = malloc((input_len + 1 + count * 5) * sizeof(char));
    if (p == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        return -1;
    }
    *output = p;

    for (i = 0; i < input_len; ++i) {
        for (esc = escapes; ; esc += 2) {
            if ((unsigned char)input[i] > **esc) {
                *(p++) = input[i];
                break;
            } else if (input[i] == **esc) {
                for (ent = *(esc + 1); *ent != '\0'; ++ent) {
                    *(p++) = *ent;
                }
                break;
            }
        }
    }

    *p = 0;
    return p - *output;
}

/*
 * Returns a copy of the given string (8-bit or Unicode) with the XML
 * control characters converted to XML character entities.
 *
 * If an 8-bit string is passed in, an 8-bit string is returned.  If a
 * Unicode string is passed in, a Unicode string is returned.
 */
static PyObject*
_escape_xml(PyObject* self, PyObject *args, const char** escapes)
{
    PyObject* input_obj;
    PyObject* input_coerce = NULL;
    PyObject* output_obj;
    char* input = NULL;
    Py_ssize_t input_len;
    char* output = NULL;
    Py_ssize_t output_len;

    if (!PyArg_ParseTuple(args, "O:escape_xml", &input_obj)) {
        return NULL;
    }

    /* First, try as Unicode */
    if (!PyBytes_Check(input_obj)) {
        input_coerce = PyObject_Str(input_obj);
    }
    if (input_coerce) {
        input = (char*)PyUnicode_AsUTF8AndSize(input_coerce, &input_len);
        if (input == NULL) {
            Py_DECREF(input_coerce);
            return NULL;
        }

        output_len = _escape_xml_impl(input, input_len, &output, escapes);
        if (output_len < 0) {
            Py_DECREF(input_coerce);
            return NULL;
        }
        if (output_len > 0) {
            Py_DECREF(input_coerce);
            output_obj = PyUnicode_FromStringAndSize(output, output_len);
            free(output);
            return output_obj;
        }
        return input_coerce;
    }

    /* Now try as bytes */
    input_coerce = PyObject_Bytes(input_obj);
    if (input_coerce) {
        if (PyBytes_AsStringAndSize(input_coerce, &input, &input_len) == -1) {
            Py_DECREF(input_coerce);
            return NULL;
        }

        output_len = _escape_xml_impl(input, input_len, &output, escapes);
        if (output_len < 0) {
            Py_DECREF(input_coerce);
            return NULL;
        }
        if (output_len > 0) {
            Py_DECREF(input_coerce);
            output_obj = PyBytes_FromStringAndSize(output, output_len);
            free(output);
            return output_obj;
        }
        return input_coerce;
    }

    PyErr_SetString(PyExc_TypeError, "must be convertible to str or bytes");
    return NULL;
}

static PyObject*
escape_xml(PyObject* self, PyObject *args)
{
    return _escape_xml(self, args, escapes);
}

static PyObject*
escape_xml_cdata(PyObject* self, PyObject *args)
{
    return _escape_xml(self, args, escapes_cdata);
}

/******************************************************************************
 * Module setup
 ******************************************************************************/

static PyMethodDef module_methods[] =
{
    {"escape_xml", (PyCFunction)escape_xml, METH_VARARGS,
     "Fast method to escape XML strings"},
    {"escape_xml_cdata", (PyCFunction)escape_xml_cdata, METH_VARARGS,
     "Fast method to escape XML strings"},
    {NULL}  /* Sentinel */
};

struct module_state {
    void* none;
};

static int module_traverse(PyObject* m, visitproc visit, void* arg)
{
    return 0;
}

static int module_clear(PyObject* m)
{
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_iterparser",
    "Fast XML parser",
    sizeof(struct module_state),
    module_methods,
    NULL,
    module_traverse,
    module_clear,
    NULL
};

PyMODINIT_FUNC
PyInit__iterparser(void)
{
    PyObject* m;
    m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    if (PyType_Ready(&IterParserType) < 0)
        return NULL;

    Py_INCREF(&IterParserType);
    PyModule_AddObject(m, "IterParser", (PyObject *)&IterParserType);

    return m;
}
