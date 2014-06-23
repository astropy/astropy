// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef TOKENIZER_H
#define TOKENIZER_H

#include <stdlib.h>

typedef enum
{
    START_LINE = 0,
    START_FIELD,
    FIELD,
    COMMENT
} tokenizer_state;

typedef enum
{
    NO_ERROR,
    INVALID_LINE,
    TOO_MANY_COLS
} err_code;

typedef struct
{
    char *source;          // single string containing all of the input
    int source_len;        // length of the input
    int source_pos;        // current index in source for tokenization
    char delimiter;        // delimiter character
    char comment;          // comment character
    char quotechar;        // quote character
    char *header_output;   // string containing header data
    char **output_cols;    // array of output strings for each column
    int output_len;        // length of each output column string
    int *row_positions;    // array of indices specifying where each row begins
    int row_pos_len;       // length of row_positions array
    int num_cols;          // number of table columns
    int num_rows;          // number of table rows
    tokenizer_state state; // current state of the tokenizer
    err_code code;         // represents the latest error that has occurred
} tokenizer_t;

/*
Example input/output
--------------------

source: "A,B,C\n10,5.,6\n1,2,3"
output_cols: ["A101", "B5.2", "C6 3"]
row_positions: [0, 1, 3]
*/

#define INITIAL_ROW_SIZE 20
#define INITIAL_COL_SIZE 50
#define INITIAL_HEADER_SIZE 50

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar);
void delete_tokenizer(tokenizer_t *tokenizer);
void resize_cols(tokenizer_t *self);
void resize_rows(tokenizer_t *self);
int tokenize(tokenizer_t *self, int line, int header);

#endif
