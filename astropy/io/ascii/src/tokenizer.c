// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include "tokenizer.h"

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar)
{
    tokenizer_t *tokenizer = (tokenizer_t *) malloc(sizeof(tokenizer_t));
    tokenizer->source = 0;
    tokenizer->source_len = 0;
    tokenizer->source_pos = 0;
    tokenizer->delimiter = delimiter;
    tokenizer->comment = comment;
    tokenizer->quotechar = quotechar;
    tokenizer->output_cols = 0;
    tokenizer->output_len = INITIAL_COL_SIZE;
    tokenizer->row_positions = 0;
    tokenizer->row_pos_len = INITIAL_ROW_SIZE;
    tokenizer->num_cols = 0;
    tokenizer->num_rows = 0;
    tokenizer->state = START_LINE;
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

void resize_cols(tokenizer_t *self)
{
}

void resize_rows(tokenizer_t *self)
{
}

// TODO: implement resize_cols, resize_rows with realloc()

#define PUSH(c) \
    if (curr_row_pos + ++field_len > self->output_len) \
	resize_cols(self); \
    self->output_cols[col][curr_row_pos + field_len - 1] = c;

#define END_FIELD() \
    if (field_len > max_field_len) \
	max_field_len = field_len; \
    field_len = 0; \
    ++col;

#define END_LINE() \
    curr_row_pos += max_field_len; \
    ++self->num_rows;
    
int tokenize(tokenizer_t *self)
{
    char c; // input character
    int col = 0; // current column
    int field_len = 0; // length of the current field
    int max_field_len = 0; // max field_len for the current row
    int curr_row_pos = 0; // first index of the current row in output_cols[0]
    self->row_positions = (int *) malloc(INITIAL_ROW_SIZE * sizeof(int));
    self->output_cols = (char **) malloc(self->num_cols * sizeof(char *));
    int i;
    for (i = 0; i < self->num_cols; ++i)
	self->output_cols[i] = (char *) calloc(1, INITIAL_COL_SIZE * sizeof(char));

    while (self->source_pos < self->source_len)
    {
	c = self->source[self->source_pos];

	switch (self->state)
	{
	case START_LINE:
	    if (c == '\n')
		break;
	    else
	    {
		if (self->num_rows >= self->row_pos_len)
		    resize_rows(self);
		self->row_positions[self->num_rows] = curr_row_pos;
		max_field_len = 0;
		col = 0;
		self->state = FIELD;
	    }
	case FIELD:
	    if (c == self->delimiter)
	    {
		END_FIELD();
	    }
	    else if (c == '\n')
	    {
		END_FIELD();
		END_LINE();
		self->state = START_LINE;
	    }
	    else
	    {
		PUSH(c);

		if (self->source_pos == self->source_len - 1)
		{
		    END_FIELD();
		    END_LINE();
		}
	    }
	    break;
	}   

	++self->source_pos;
    }

    return 0;
}
