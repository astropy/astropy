// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include "tokenizer.h"

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar)
{
    tokenizer_t *tokenizer = (tokenizer_t *) calloc(1, sizeof(tokenizer_t));
    tokenizer->source = 0;
    tokenizer->source_len = 0;
    tokenizer->source_pos = 0;
    tokenizer->delimiter = delimiter;
    tokenizer->comment = comment;
    tokenizer->quotechar = quotechar;
    tokenizer->output_cols = 0;
    tokenizer->row_positions = 0;
    tokenizer->num_cols = 0;
    tokenizer->err_code = 0; //TODO: make const int ERR_* = ... declarations for clarity
    return tokenizer;
}

void delete_tokenizer(tokenizer_t *tokenizer)
{
    // TODO: call free() on tokenizer elements
    free(tokenizer);
}

int tokenize_header(tokenizer_t *self)
{
    return 0;
}

int tokenize(tokenizer_t *self)
{
    return 0;
}
