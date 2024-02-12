// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/types.h>

#ifdef _MSC_VER
    #define inline __inline
    #ifndef NAN
        static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
        #define NAN (*(const double *) __nan)
    #endif
    #ifndef INFINITY
        static const unsigned long __infinity[2] = {0x00000000, 0x7ff00000};
        #define INFINITY (*(const double *) __infinity)
    #endif
#else
    #ifndef INFINITY
        #define INFINITY (1.0/0.0)
    #endif
    #ifndef NAN
        #define NAN (INFINITY-INFINITY)
    #endif
#endif

typedef enum
{
    START_LINE = 0,
    START_FIELD,
    START_QUOTED_FIELD,
    FIELD,
    QUOTED_FIELD,
    QUOTED_FIELD_DOUBLE_QUOTE,
    COMMENT,
} tokenizer_state;

typedef enum
{
    NO_ERROR,
    INVALID_LINE,
    TOO_MANY_COLS,
    NOT_ENOUGH_COLS,
    CONVERSION_ERROR,
    OVERFLOW_ERROR
} err_code;

typedef struct
{
    char *source;          // single string containing all of the input
    size_t source_len;      // length of the input
    size_t source_pos;      // current index in source for tokenization
    char delimiter;        // delimiter character
    char comment;          // comment character
    char quotechar;        // quote character
    char expchar;          // exponential character in scientific notation
    char newline;          // EOL character
    char **output_cols;    // array of output strings for each column
    char **col_ptrs;       // array of pointers to current output position for each col
    size_t *output_len;    // length of each output column string
    int num_cols;          // number of table columns
    int num_rows;          // number of table rows
    int fill_extra_cols;   // represents whether or not to fill rows with too few values
    tokenizer_state state; // current state of the tokenizer
    err_code code;         // represents the latest error that has occurred
    int iter_col;          // index of the column being iterated over
    char *curr_pos;        // current iteration position
    char *buf;             // buffer for empty data
    int strip_whitespace_lines;  // whether to strip whitespace at the beginning and end of lines
    int strip_whitespace_fields; // whether to strip whitespace at the beginning and end of fields
    int use_fast_converter;      // whether to use the fast converter for floats
    char *comment_lines;   // single null-delimited string containing comment lines
    int comment_lines_len; // length of comment_lines in memory
    int comment_pos;       // current index in comment_lines
} tokenizer_t;

/*
Example input/output
--------------------

source: "A,B,C\n10,5.,6\n1,2,3"
output_cols: ["A\x0010\x001", "B\x005.\x002", "C\x006\x003"]
*/

#define INITIAL_COL_SIZE 500
#define INITIAL_COMMENT_LEN 50

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar, char expchar,
                              int fill_extra_cols, int strip_whitespace_lines,
                              int strip_whitespace_fields, int use_fast_converter);
void delete_tokenizer(tokenizer_t *tokenizer);
void delete_data(tokenizer_t *tokenizer);
void resize_col(tokenizer_t *self, int index);
void resize_comments(tokenizer_t *self);
int skip_lines(tokenizer_t *self, int offset, int header);
int tokenize(tokenizer_t *self, int end, int header, int num_cols);
int64_t str_to_int64_t(tokenizer_t *self, char *str);
double str_to_double(tokenizer_t *self, char *str);
double xstrtod(const char *str, char **endptr, char decimal,
               char expchar, char tsep, int skip_trailing);
void start_iteration(tokenizer_t *self, int col);
char *next_field(tokenizer_t *self, int *size);
long file_len(FILE *fhandle);
char *get_line(char *ptr, size_t *len, size_t map_len);
void reset_comments(tokenizer_t *self);

#endif
