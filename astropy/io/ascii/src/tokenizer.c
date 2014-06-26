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
    tokenizer->col_ptrs = 0;
    tokenizer->header_len = 0;
    tokenizer->output_len = 0;
    tokenizer->num_cols = 0;
    tokenizer->num_rows = 0;
    tokenizer->state = START_LINE;
    tokenizer->code = NO_ERROR;
    return tokenizer;
}

void delete_data(tokenizer_t *tokenizer)
{
    // Don't free tokenizer->source because it points to part of an already freed Python object
    free(tokenizer->header_output);
    int i;
    if (tokenizer->output_cols)
	for (i = 0; i < tokenizer->num_cols; ++i)
	    free(tokenizer->output_cols[i]);

    free(tokenizer->output_cols);
    free(tokenizer->col_ptrs);
    tokenizer->header_output = 0;
    tokenizer->output_cols = 0;
    tokenizer->col_ptrs = 0;
}

void delete_tokenizer(tokenizer_t *tokenizer)
{
    delete_data(tokenizer);
    free(tokenizer);
}

void resize_col(tokenizer_t *self, int index)
{
    self->output_len[index] *= 2;
    self->output_cols[index] = (char *) realloc(self->output_cols[index], self->output_len[index] * sizeof(char));
}

#define PUSH(c)								\
    if (header)								\
    {									\
        if (output_pos >= self->header_len)				\
        {								\
            self->header_len *= 2;					\
            self->header_output = (char *) realloc(self->header_output, \
				      self->header_len * sizeof(char)); \
        }								\
	self->header_output[output_pos++] = c;				\
    }									\
    else if (col < self->num_cols && use_cols[real_col])		\
    {									\
	if (self->col_ptrs[col] - self->output_cols[col] >= self->output_len[col]) \
	{								\
	    resize_col(self, col);					\
	}								\
        *self->col_ptrs[col]++ = c;					\
    }

#define END_FIELD()							\
    if (header || use_cols[real_col])					\
    {									\
	PUSH('\x00');							\
	if (!header && ++col > self->num_cols)				\
	    RETURN(TOO_MANY_COLS);					\
    }									\
    ++real_col;

#define END_LINE()				\
    if (header)					\
	done = 1;				\
    else					\
    {						\
	while (col < self->num_cols)		\
	{					\
            PUSH(' ');				\
	    END_FIELD();			\
	}					\
    }						\
    ++self->num_rows;				\

#define RETURN(c) { self->code = c; return c; }

int tokenize(tokenizer_t *self, int line, int header, int *use_cols)
{
    delete_data(self); // clear old reading data
    char c; // input character
    int col = 0; // current column ignoring possibly excluded columns
    int real_col = 0; // current column taking excluded columns into account
    int output_pos = 0; // current position in header output string
    self->header_len = INITIAL_HEADER_SIZE;
    self->source_pos = 0;
    self->num_rows = 0;
    int i = 0;
    int empty = 1;
    int comment = 0;

    //TODO: different error for no data
    //TODO: decide what to do about whitespace delimiter here
    while (i < line)
    {
	if (self->source_pos >= self->source_len - 1) // ignore final newline
	    RETURN(INVALID_LINE);
	if (self->source[self->source_pos] != '\n' && empty)
	{
	    empty = 0;
	    if (self->source[self->source_pos] == self->comment)
		comment = 1;
	}
	else if (self->source[self->source_pos++] == '\n')
	{
	    if (!empty && !comment)
		++i;
	    empty = 1;
	    comment = 0;
	}
    }

    if (header)
	self->header_output = (char *) calloc(1, INITIAL_HEADER_SIZE * sizeof(char));
    else
    {
	self->output_cols = (char **) malloc(self->num_cols * sizeof(char *));
	self->col_ptrs = (char **) malloc(self->num_cols * sizeof(char *));
	self->output_len = (int *) malloc(self->num_cols * sizeof(int));

	for (i = 0; i < self->num_cols; ++i)
	{
	    self->output_cols[i] = (char *) calloc(1, INITIAL_COL_SIZE * sizeof(char));
	    self->col_ptrs[i] = self->output_cols[i];
	    self->output_len[i] = INITIAL_COL_SIZE;
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
	    if (c == '\n' || c == ' ' || c == '\t') // TODO: make an option not to strip whitespace (for tab-delimited, etc.)
		break;
	    else if (c == self->comment)
	    {
		self->state = COMMENT;
		break;
	    }
	    col = 0;
	    real_col = 0;
	    self->state = START_FIELD;
	
	case START_FIELD:
	    if (c == ' ' || c == '\t') // TODO: strip whitespace at the end of fields as well
		break;
	    else if (c == self->delimiter)
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
	    else
	    {
		PUSH(c);
	    }
	    break;

	case COMMENT:
	    if (c == '\n')
		self->state = START_LINE;
	    break;
	}   

	++self->source_pos;
    }

    RETURN(0);
}

int int_size()
{
    return 8 * sizeof(int);
}
