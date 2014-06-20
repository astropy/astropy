// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include "tokenizer.h"

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar)
{
    tokenizer_t *tokenizer = (tokenizer_t *) calloc(1, sizeof(tokenizer_t));
    tokenizer->delimiter = delimiter;
    tokenizer->comment = comment;
    tokenizer->quotechar = quotechar;
}

int tokenize(tokenizer_t *self, int rows)
{
    
