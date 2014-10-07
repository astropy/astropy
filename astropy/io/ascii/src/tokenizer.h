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

typedef enum
{
    START_LINE = 0,
    START_FIELD,
    START_QUOTED_FIELD,
    FIELD,
    QUOTED_FIELD,
    QUOTED_FIELD_NEWLINE,
    COMMENT,
    CARRIAGE_RETURN
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
    int source_len;        // length of the input
    int source_pos;        // current index in source for tokenization
    char delimiter;        // delimiter character
    char comment;          // comment character
    char quotechar;        // quote character
    char **output_cols;    // array of output strings for each column
    char **col_ptrs;       // array of pointers to current output position for each col
    int *output_len;       // length of each output column string
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
} tokenizer_t;

/*
Example input/output
--------------------

source: "A,B,C\n10,5.,6\n1,2,3"
output_cols: ["A\x0010\x001", "B\x005.\x002", "C\x006\x003"]
*/

#define INITIAL_COL_SIZE 500

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar, int fill_extra_cols,
                              int strip_whitespace_lines, int strip_whitespace_fields,
                              int use_fast_converter);
void delete_tokenizer(tokenizer_t *tokenizer);
void delete_data(tokenizer_t *tokenizer);
void resize_col(tokenizer_t *self, int index);
int skip_lines(tokenizer_t *self, int offset, int header);
int tokenize(tokenizer_t *self, int end, int header, int num_cols);
long str_to_long(tokenizer_t *self, char *str);
double str_to_double(tokenizer_t *self, char *str);
double xstrtod(const char *str, char **endptr, char decimal,
               char sci, char tsep, int skip_trailing);
void start_iteration(tokenizer_t *self, int col);
int finished_iteration(tokenizer_t *self);
char *next_field(tokenizer_t *self, int *size);
long file_len(FILE *fhandle);
char *get_line(char *ptr, int *len, int map_len);

#endif
