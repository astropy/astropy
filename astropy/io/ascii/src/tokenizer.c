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
    tokenizer->header_output = 0;
    tokenizer->output_cols = 0;
    tokenizer->output_len = INITIAL_COL_SIZE;
    tokenizer->row_positions = 0;
    tokenizer->row_pos_len = INITIAL_ROW_SIZE;
    tokenizer->num_cols = 0;
    tokenizer->num_rows = 0;
    tokenizer->state = START_LINE;
    tokenizer->code = NO_ERROR;
    return tokenizer;
}

void delete_tokenizer(tokenizer_t *tokenizer)
{
    // Don't free tokenizer->source because it points to part of an already freed Python object
    free(tokenizer->row_positions);
    free(tokenizer->header_output);
    int i;
    for (i = 0; i < tokenizer->num_cols; ++i)
	free(tokenizer->output_cols[i]);
    free(tokenizer->output_cols);
    free(tokenizer);
}

void resize_cols(tokenizer_t *self)
{
    int i;
    self->output_len *= 2;
    for (i = 0; i < self->num_cols; ++i)
	self->output_cols[i] = (char *)realloc(self->output_cols[i], self->output_len * sizeof(char));
}

void resize_rows(tokenizer_t *self)
{
    self->row_pos_len *= 2;
    self->row_positions = (int *)realloc(self->row_positions, self->row_pos_len * sizeof(int));
}

#define PUSH(c) \
    if (header) \
    { \
        if (output_pos >= output_len) \
        { \
            output_len *= 2; \
            self->header_output = (char *) realloc(self->header_output, output_len * sizeof(char)); \
        } \
	self->header_output[output_pos++] = c; \
    } \
    else \
    { \
        if (curr_row_pos + ++field_len > self->output_len) \
            resize_cols(self); \
        self->output_cols[col][curr_row_pos + field_len - 1] = c; \
    }

#define END_FIELD() \
    if (header) \
    { \
	PUSH('\x00'); \
    } \
    else \
    { \
        if (field_len > max_field_len) \
	    max_field_len = field_len; \
        field_len = 0; \
        if (++col > self->num_cols) \
	    RETURN(TOO_MANY_COLS); \
    }

#define END_LINE() \
    if (header) \
	done = 1; \
    else \
    { \
	while (col < self->num_cols) \
	{ \
            PUSH(' '); \
	    END_FIELD(); \
	} \
    } \
    curr_row_pos += max_field_len; \
    ++self->num_rows; \

#define RETURN(c) { self->code = c; return c; }

int tokenize(tokenizer_t *self, int line, int header)
{
    char c; // input character
    int col = 0; // current column
    int field_len = 0; // length of the current field
    int max_field_len = 0; // max field_len for the current row
    int curr_row_pos = 0; // first index of the current row in output_cols[0]
    int output_len = INITIAL_HEADER_SIZE;
    int output_pos = 0;
    int i = 0;
    self->source_pos = 0;
    self->num_rows = 0;

    //TODO: fix this for empty/comment lines, different error for no data, etc.
    while (i < line)
    {
	if (self->source_pos >= self->source_len)
	    RETURN(INVALID_LINE);
	if (self->source[self->source_pos++] == '\n')
	    ++i;
    }

    if (header)
	self->header_output = (char *) calloc(1, INITIAL_HEADER_SIZE * sizeof(char));
    else
    {
	if (!self->row_positions)
	    self->row_positions = (int *) malloc(INITIAL_ROW_SIZE * sizeof(int));
	if (!self->output_cols)
	{
	    self->output_cols = (char **) malloc(self->num_cols * sizeof(char *));
	    int i;
	    for (i = 0; i < self->num_cols; ++i)
		self->output_cols[i] = (char *) calloc(1, INITIAL_COL_SIZE * sizeof(char));
	}
    }

    int done = 0;
    self->state = START_LINE;

    while (self->source_pos < self->source_len && !done)
    {
	c = self->source[self->source_pos];

	switch (self->state)
	{
	case START_LINE:
	    if (c == '\n')
		break;
	    else if (c == self->comment)
	    {
		self->state = COMMENT_LINE;
		break;
	    }
	    else if (!header)
	    {
		if (self->num_rows >= self->row_pos_len)
		    resize_rows(self);
		self->row_positions[self->num_rows] = curr_row_pos;
		max_field_len = 0;
		col = 0;
	    }
	    self->state = START_FIELD;
	
	case START_FIELD:
	    if (c == self->delimiter)
	    {
		PUSH(' ');
		END_FIELD();
		break;
	    }
	    else if (c == '\n')
	    {
		END_LINE();
		self->state = START_LINE;
		break;
	    }
	    else if (c == self->comment)
	    {
		self->state = COMMENT;
		break;
	    }
	    self->state = FIELD;

	case FIELD:
	    if (c == self->delimiter)
	    {
		END_FIELD();
		self->state = START_FIELD;
	    }
	    else if (c == '\n')
	    {
		END_FIELD();
		END_LINE();
		self->state = START_LINE;
	    }
	    else if (c == self->comment)
	    {
		END_FIELD();
		self->state = COMMENT;
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

	case COMMENT: //TODO: figure out whether this should be here (old readers don't support non-line comments)
	    if (c == '\n')
	    {
		END_LINE();
		self->state = START_LINE;
	    }
	    break;

	case COMMENT_LINE:
	    if (c == '\n')
		self->state = START_LINE;
	    break;
	}   

	++self->source_pos;
    }

    RETURN(0);
}
