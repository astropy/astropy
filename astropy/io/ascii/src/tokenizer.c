// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include "tokenizer.h"

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar, int fill_extra_cols,
                              int strip_whitespace_lines, int strip_whitespace_fields)
{
    // Create the tokenizer in memory
    tokenizer_t *tokenizer = (tokenizer_t *) malloc(sizeof(tokenizer_t));

    // Initialize the tokenizer fields
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
    tokenizer->fill_extra_cols = fill_extra_cols;
    tokenizer->state = START_LINE;
    tokenizer->code = NO_ERROR;
    tokenizer->iter_col = 0;
    tokenizer->curr_pos = 0;
    tokenizer->strip_whitespace_lines = strip_whitespace_lines;
    tokenizer->strip_whitespace_fields = strip_whitespace_fields;

    // This is a bit of a hack -- buf holds an empty string to represent empty field values
    tokenizer->buf = calloc(2, sizeof(char));

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

    // Set pointers to 0 so we don't use freed memory when reading over again
    tokenizer->header_output = 0;
    tokenizer->output_cols = 0;
    tokenizer->col_ptrs = 0;
}

void delete_tokenizer(tokenizer_t *tokenizer)
{
    delete_data(tokenizer);
    free(tokenizer->buf);
    free(tokenizer);
}

void resize_col(tokenizer_t *self, int index)
{
    // Temporarily store the position in output_cols[index] to which col_ptrs[index] points
    int diff = self->col_ptrs[index] - self->output_cols[index];

    // Double the size of the column string
    self->output_cols[index] = (char *) realloc(self->output_cols[index], 2 * self->output_len[index] * sizeof(char));

    // Set the second (newly allocated) half of the column string to all zeros
    memset(self->output_cols[index] + self->output_len[index] * sizeof(char), 0, self->output_len[index] * sizeof(char));

    self->output_len[index] *= 2;
    // realloc() might move the address in memory, so we have to move col_ptrs[index] to an offset of the new address
    self->col_ptrs[index] = self->output_cols[index] + diff;
}

void resize_header(tokenizer_t *self)
{
    // Double the size of the header string
    self->header_output = (char *) realloc(self->header_output, 2 * self->header_len * sizeof(char));

    // Set the second half of the header string to all zeros
    memset(self->header_output + self->header_len * sizeof(char), 0, self->header_len * sizeof(char));
    self->header_len *= 2;
}

/*
  If we are tokenizing the header, resize the header string if necessary 
  and then set the char at the current output position to c, incrementing
  the current position afterwards.

  If we are tokenizing data but we have too many columns or use_cols[real_col]
  is 0 (meaning that this column should be excluded from output), do nothing.
  Otherwise, resize the column string if necessary and then append c to the
  end of the column string, incrementing the column position pointer.
*/

#define PUSH(c)								\
    if (header)								\
    {									\
        if (output_pos >= self->header_len)				\
        {								\
	    resize_header(self);					\
	}								\
	self->header_output[output_pos++] = c;				\
    }									\
    else if (col < self->num_cols && use_cols[real_col])		\
    {									\
	if (self->col_ptrs[col] - self->output_cols[col] >= self->output_len[col] * sizeof(char)) \
	{								\
	    resize_col(self, col);					\
	}								\
        *self->col_ptrs[col]++ = c;					\
    }

/*
  Unless this column will be excluded from output, append a null
  byte to the end of the column string as a field delimiting marker.
  Increment the variable col and eturn the value TOO_MANY_COLS if
  there are too many columns in this row. Increment real_col even if
  this column should be ignored.
*/

#define END_FIELD()							\
    if (header || real_col >= use_cols_len || use_cols[real_col])	\
    {									\
	PUSH('\x00');							\
	if (!header && ++col > self->num_cols)				\
	    RETURN(TOO_MANY_COLS);					\
    }									\
    ++real_col;

/*
  If we are tokenizing the header, end after the first line.
  Handle the possibility of insufficient columns appropriately;
  if fill_extra_cols=1, then append empty fields, but otherwise
  return an error. Increment our row count and possibly end if
  all the necessary rows have already been parsed.
*/

#define END_LINE()					\
    if (header)						\
	done = 1;					\
    else if (self->fill_extra_cols)			\
	while (col < self->num_cols)			\
	{						\
            PUSH('\x01');				\
	    END_FIELD();				\
	}						\
    else if (col < self->num_cols)			\
	RETURN(NOT_ENOUGH_COLS);			\
    ++self->num_rows;					\
    if (end != -1 && self->num_rows == end - start)	\
	done = 1;

// Set the error code to c for later retrieval and return c
#define RETURN(c) { self->code = c; return c; }

int tokenize(tokenizer_t *self, int start, int end, int header, int *use_cols, int use_cols_len)
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
    
    //TODO: fix quoting issues here...maybe we can tokenize these rows in a new non-output mode
    while (i < start)
    {
	if (self->source_pos >= self->source_len - 1) // ignore final newline
	    RETURN(INVALID_LINE);
	if (self->source[self->source_pos] != '\n' && empty) // first significant character encountered
	{
	    empty = 0;
	    if (self->source[self->source_pos] == self->comment) // comment line
		comment = 1;
	}
	else if (self->source[self->source_pos++] == '\n')
	{
	    if (!empty && !comment) // significant line
		++i;
            // Start by assuming a line is empty and non-commented
	    empty = 1;
	    comment = 0;
	}
    }
    
    // Allocate memory for structures used during tokenization
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
            // Make each col_ptrs pointer point to the beginning of the column string
	    self->col_ptrs[i] = self->output_cols[i];
	    self->output_len[i] = INITIAL_COL_SIZE;
	}
    }
    
    int done = (end != -1 && end <= start); // Make sure the parameter end is valid
    int repeat;
    self->state = START_LINE;
    
    // Loop until all of source has been read or we finish for some other reason
    while (self->source_pos < self->source_len && !done)
    {
	c = self->source[self->source_pos]; // next character of input
	repeat = 1;
	
        // Keep parsing the same character
	while (repeat && !done)
	{
	    repeat = 0; // don't repeat by default; handling might change this

	    switch (self->state)
	    {
	    case START_LINE:
                if (c == '\n')
                    break;
                else if ((c == ' ' || c == '\t') && self->strip_whitespace_lines)
                    break;                    
		else if (c == self->comment) // comment line; ignore
		{
		    self->state = COMMENT;
		    break;
		}
                // initialize variables for the beginning of line parsing
		col = 0;
		real_col = 0;
		self->state = START_FIELD;
		repeat = 1; // parse the same character again in mode START_FIELD
		break;
		
	    case START_FIELD:
		if ((c == ' ' || c == '\t') && self->strip_whitespace_fields) // strip whitespace before field begins
		    ;
		else if (c == self->delimiter) // field ends before it begins
		{
		    PUSH('\x01'); // indicates empty field
		    END_FIELD();
		}
		else if (c == '\n') // TODO: right stripping of fields
		{
                    if (self->strip_whitespace_lines)
                    {
                        // Move on if the delimiter is whitespace, e.g. '1 2 3   '->['1','2','3']
                        if (self->delimiter == ' ' || self->delimiter == '\t')
                        {
                            END_FIELD();
                        }
                        // Register an empty field if non-whitespace delimiter, e.g. '1,2, '->['1','2','']
                        else
                        {
                            PUSH('\x01');
                            END_FIELD();
                        }                            
                    }

                    else if (!self->strip_whitespace_lines)
                    {
                        // In this case we don't want to left-strip the field, so we backtrack
                        int tmp = self->source_pos;
                        --self->source_pos;

                        while (self->source_pos >= 0 && self->source[self->source_pos] != self->delimiter
                               && self->source[self->source_pos] != '\n')
                        {
                            --self->source_pos;
                        }

                        // backtracked to line beginning
                        if (self->source_pos == -1 || self->source[self->source_pos] == '\n')
                            self->source_pos = tmp;
                        else
                        {
                            ++self->source_pos;

                            if (self->source_pos == tmp)
                            {
                                PUSH('\x01'); // no whitespace, just an empty field
                            }

                            else
                                while (self->source_pos < tmp)
                                {
                                    PUSH(self->source[self->source_pos]); // append whitespace characters
                                    ++self->source_pos;
                                }

                            END_FIELD(); // whitespace counts as a field
                        }
                    }
                        
		    END_LINE();
		    self->state = START_LINE;
		}
		else if (c == self->quotechar) // start parsing quoted field
		    self->state = START_QUOTED_FIELD;
		else
		{
                    // Valid field character, parse again in FIELD mode
		    repeat = 1;
		    self->state = FIELD;
		}
		break;
		
	    case START_QUOTED_FIELD:
		if (c == ' ' || c == '\t') // ignore initial whitespace
		    ;
		else if (c == self->quotechar) // empty quotes
		{
		    PUSH('\x01'); // indicates empty field
		    END_FIELD();
		}
		else
		{
                    // Valid field character, parse again in QUOTED_FIELD mode
		    self->state = QUOTED_FIELD;
		    repeat = 1;
		}
		break;
		
	    case FIELD:
		if (c == self->delimiter)
		{
                    // End of field, look for new field
		    END_FIELD();
		    self->state = START_FIELD;
		}
		else if (c == '\n')
		{
                    // Line ending, stop parsing both field and line
		    END_FIELD();
		    END_LINE();
		    self->state = START_LINE;
		}
		else
		{                    
		    PUSH(c);
		}
		break;
		
	    case QUOTED_FIELD:
		if (c == self->quotechar) // Parse rest of field normally, e.g. "ab"c
		    self->state = FIELD;
		else if (c == '\n')
		    self->state = QUOTED_FIELD_NEWLINE;
		else
		{
		    PUSH(c);
		}
		break;

	    case QUOTED_FIELD_NEWLINE:
                // Ignore initial whitespace/newlines
		if (c == ' ' || c == '\t' || c == '\n')
		    ;
		else if (c == self->quotechar)
		    self->state = FIELD; // TODO: fix this for empty data
		else
		{
                    // Once data begins, parse it as a normal quoted field
		    repeat = 1;
		    self->state = QUOTED_FIELD;
		}
		break;
		
	    case COMMENT:
		if (c == '\n')
		    self->state = START_LINE;
		break; // keep looping until we find a newline
	    }
	}
	
	++self->source_pos;
    }
    
    RETURN(0);
}

int int_size(void)
{
    return 8 * sizeof(int);
}

long str_to_long(tokenizer_t *self, char *str)
{
    errno = 0;
    char *tmp;
    long ret = strtol(str, &tmp, 0);

    if (tmp == str || *tmp != '\0')
 	self->code = CONVERSION_ERROR;
    else if (errno == ERANGE)
        self->code = OVERFLOW_ERROR;

    return ret;
}

double str_to_double(tokenizer_t *self, char *str)
{
    char *tmp;
    double ret = strtod(str, &tmp);

    if (tmp == str || *tmp != 0) // TODO: make sure this is right
	self->code = CONVERSION_ERROR;
    else if (errno == ERANGE)
        self->code = OVERFLOW_ERROR;

    return ret;
}

void start_iteration(tokenizer_t *self, int col)
{
    // Begin looping over the column string with index col
    self->iter_col = col;
    self->curr_pos = self->output_cols[col]; // Start at the initial pointer position
}

int finished_iteration(tokenizer_t *self)
{
    // Iteration is finished if we've exceeded the length of the column string or the pointer is at an empty byte
    return (self->curr_pos - self->output_cols[self->iter_col] >= self->output_len[self->iter_col] * sizeof(char)
	    || *self->curr_pos == '\x00');
}

char *next_field(tokenizer_t *self)
{
    char *tmp = self->curr_pos;

    while (*self->curr_pos != '\x00') // pass through the entire field until reaching the delimiter
	++self->curr_pos;

    ++self->curr_pos; // next field begins after the delimiter
    if (*tmp == '\x01') // empty field; this is a hack
	return self->buf;
    else
	return tmp;
}
