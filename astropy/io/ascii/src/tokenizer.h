// Licensed under a 3-clause BSD style license - see LICENSE.rst

#ifndef TOKENIZER_H
#define TOKENIZER_H

typedef struct
{
    char *source;          // single string containing all of the input
    int source_len;        // length of the input
    char delimiter;        // delimiter character
    char comment;          // comment character
    char quotechar;        // quote character
    char **output_cols;    // array of output strings for each column
    int *row_positions;    // array of indices specifying where each row begins
} tokenizer_t;

/ *
Example input/output
--------------------

source: "A,B,C\n10,5.,6\n1,2,3"
output_cols: ["A101", "B5.2", "C6 3"]
row_positions: [0, 1, 3]
* /

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar);
int tokenize(tokenizer_t *self, int rows);

#endif
