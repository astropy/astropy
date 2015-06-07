// Licensed under a 3-clause BSD style license - see LICENSE.rst

#include "tokenizer.h"

tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar,
                              int fill_extra_cols, int strip_whitespace_lines,
                              int strip_whitespace_fields, int use_fast_converter)
{
    // Create the tokenizer in memory
    tokenizer_t *tokenizer = (tokenizer_t *) malloc(sizeof(tokenizer_t));

    // Initialize the tokenizer fields
    tokenizer->source = NULL;
    tokenizer->source_len = 0;
    tokenizer->source_pos = 0;
    tokenizer->delimiter = delimiter;
    tokenizer->comment = comment;
    tokenizer->quotechar = quotechar;
    tokenizer->output_cols = NULL;
    tokenizer->col_ptrs = NULL;
    tokenizer->output_len = NULL;
    tokenizer->num_cols = 0;
    tokenizer->num_rows = 0;
    tokenizer->fill_extra_cols = fill_extra_cols;
    tokenizer->state = START_LINE;
    tokenizer->code = NO_ERROR;
    tokenizer->iter_col = 0;
    tokenizer->curr_pos = NULL;
    tokenizer->strip_whitespace_lines = strip_whitespace_lines;
    tokenizer->strip_whitespace_fields = strip_whitespace_fields;
    tokenizer->use_fast_converter = use_fast_converter;
    tokenizer->comment_lines = (char *) malloc(INITIAL_COMMENT_LEN);
    tokenizer->comment_pos = 0;
    tokenizer->comment_lines_len = 0;

    // This is a bit of a hack -- buf holds an empty string to represent
    // empty field values
    tokenizer->buf = calloc(2, sizeof(char));

    return tokenizer;
}


void delete_data(tokenizer_t *tokenizer)
{
    // Don't free tokenizer->source because it points to part of
    // an already freed Python object
    int i;

    if (tokenizer->output_cols)
    {
        for (i = 0; i < tokenizer->num_cols; ++i)
        {
            free(tokenizer->output_cols[i]);
        }
    }

    free(tokenizer->output_cols);
    free(tokenizer->col_ptrs);
    free(tokenizer->output_len);

    // Set pointers to 0 so we don't use freed memory when reading over again
    tokenizer->output_cols = 0;
    tokenizer->col_ptrs = 0;
    tokenizer->output_len = 0;
}


void delete_tokenizer(tokenizer_t *tokenizer)
{
    delete_data(tokenizer);
    free(tokenizer->comment_lines);
    free(tokenizer->buf);
    free(tokenizer);
}


void resize_col(tokenizer_t *self, int index)
{
    // Temporarily store the position in output_cols[index] to
    // which col_ptrs[index] points
    long diff = self->col_ptrs[index] - self->output_cols[index];

    // Double the size of the column string
    self->output_cols[index] = (char *) realloc(self->output_cols[index], 2 *
                                                self->output_len[index] * sizeof(char));

    // Set the second (newly allocated) half of the column string to all zeros
    memset(self->output_cols[index] + self->output_len[index] * sizeof(char), 0,
           self->output_len[index] * sizeof(char));

    self->output_len[index] *= 2;
    // realloc() might move the address in memory, so we have to move
    // col_ptrs[index] to an offset of the new address
    self->col_ptrs[index] = self->output_cols[index] + diff;
}


void resize_comments(tokenizer_t *self)
{
    // Double the size of the comments string
    self->comment_lines = (char *) realloc(self->comment_lines,
                                           self->comment_pos + 1);
    // Set the second (newly allocated) half of the column string to all zeros
    memset(self->comment_lines + self->comment_lines_len * sizeof(char), 0,
           (self->comment_pos + 1 - self->comment_lines_len) * sizeof(char));

    self->comment_lines_len = self->comment_pos + 1;
}

/*
  Resize the column string if necessary and then append c to the
  end of the column string, incrementing the column position pointer.
*/
static inline void push(tokenizer_t *self, char c, int col)
{
    if (self->col_ptrs[col] - self->output_cols[col] >=
        self->output_len[col])
    {
        resize_col(self, col);
    }

    *self->col_ptrs[col]++ = c;
}


/*
  Resize the comment string if necessary and then append c to the
  end of the comment string.
*/
static inline void push_comment(tokenizer_t *self, char c)
{
    if (self->comment_pos >= self->comment_lines_len)
    {
        resize_comments(self);
    }
    self->comment_lines[self->comment_pos++] = c;
}


static inline void end_comment(tokenizer_t *self)
{
    // Signal empty comment by inserting \x01
    if (self->comment_pos == 0 || self->comment_lines[self->comment_pos - 1] == '\x00')
    {
        push_comment(self, '\x01');
    }
    push_comment(self, '\x00');
}


#define PUSH(c) push(self, c, col)


/* Set the state to START_FIELD and begin with the assumption that
   the field is entirely whitespace in order to handle the possibility
   that the comment character is found before any non-whitespace even
   if whitespace stripping is disabled.
*/
#define BEGIN_FIELD()                           \
    self->state = START_FIELD;                  \
    whitespace = 1


/*
  First, backtrack to eliminate trailing whitespace if strip_whitespace_fields
  is true. If the field is empty, push '\x01' as a marker.
  Append a null byte to the end of the column string as a field delimiting marker.
  Increment the variable col if we are tokenizing data.
*/
static inline void end_field(tokenizer_t *self, int *col, int header)
{
    if (self->strip_whitespace_fields &&
            self->col_ptrs[*col] != self->output_cols[*col])
    {
        --self->col_ptrs[*col];
        while (*self->col_ptrs[*col] == ' ' || *self->col_ptrs[*col] == '\t')
        {
            *self->col_ptrs[*col]-- = '\x00';
        }
        ++self->col_ptrs[*col];
    }
    if (self->col_ptrs[*col] == self->output_cols[*col] ||
            self->col_ptrs[*col][-1] == '\x00')
    {
        push(self, '\x01', *col);
    }
    push(self, '\x00', *col);
    if (!header) {
        ++*col;
    }
}


#define END_FIELD() end_field(self, &col, header)


// Set the error code to c for later retrieval and return c
#define RETURN(c)                                               \
    do {                                                        \
        self->code = c;                                         \
        return c;                                               \
    } while (0)


/*
  If we are tokenizing the header, end after the first line.
  Handle the possibility of insufficient columns appropriately;
  if fill_extra_cols=1, then append empty fields, but otherwise
  return an error. Increment our row count and possibly end if
  all the necessary rows have already been parsed.
*/
static inline int end_line(tokenizer_t *self, int col, int header, int end,
                           tokenizer_state *old_state)
{
    if (header)
    {
        ++self->source_pos;
        RETURN(NO_ERROR);
    }
    else if (self->fill_extra_cols)
    {
        while (col < self->num_cols)
        {
                PUSH('\x01');
            END_FIELD();
        }
    }
    else if (col < self->num_cols)
    {
        RETURN(NOT_ENOUGH_COLS);
    }

    ++self->num_rows;
    *old_state = START_LINE;

    if (end != -1 && self->num_rows == end)
    {
        ++self->source_pos;
        RETURN(NO_ERROR);
    }
    return -1;
}


#define END_LINE() if (end_line(self, col, header, end, &old_state) != -1) return self->code


#define HANDLE_CR() old_state = self->state; self->state = CARRIAGE_RETURN


int skip_lines(tokenizer_t *self, int offset, int header)
{
    int signif_chars = 0;
    int comment = 0;
    int i = 0;
    char c;

    while (i < offset)
    {
        if (self->source_pos >= self->source_len)
        {
            if (header)
                RETURN(INVALID_LINE); // header line is required
            else
                RETURN(NO_ERROR); // no data in input
        }

        c = self->source[self->source_pos];

        if (c == '\r' || c == '\n')
        {
            if (c == '\r' && self->source_pos < self->source_len - 1 &&
                self->source[self->source_pos + 1] == '\n')
            {
                ++self->source_pos; // skip \n in \r\n
            }
            if (!comment && signif_chars > 0)
                ++i;
            else if (comment && !header)
                end_comment(self);
            // Start by assuming a line is empty and non-commented
            signif_chars = 0;
            comment = 0;
        }
        else if ((c != ' ' && c != '\t') || !self->strip_whitespace_lines)
        {
                // comment line
                if (!signif_chars && self->comment != 0 && c == self->comment)
                    comment = 1;
                else if (comment && !header)
                    push_comment(self, c);

                // significant character encountered
                ++signif_chars;
        }
        else if (comment && !header)
        {
            push_comment(self, c);
        }

            ++self->source_pos;
    }

    RETURN(NO_ERROR);
}


int tokenize(tokenizer_t *self, int end, int header, int num_cols)
{
    char c; // input character
    int col = 0; // current column ignoring possibly excluded columns
    tokenizer_state old_state = START_LINE; // last state the tokenizer was in before CR mode
    int parse_newline = 0; // explicit flag to treat current char as a newline
    int i = 0;
    int whitespace = 1;
    delete_data(self); // clear old reading data
    self->num_rows = 0;
    self->comment_lines_len = INITIAL_COMMENT_LEN;

    if (header)
        self->num_cols = 1; // store header output in one column
    else
        self->num_cols = num_cols;

    // Allocate memory for structures used during tokenization
    self->output_cols = (char **) malloc(self->num_cols * sizeof(char *));
    self->col_ptrs = (char **) malloc(self->num_cols * sizeof(char *));
    self->output_len = (int *) malloc(self->num_cols * sizeof(int));

    for (i = 0; i < self->num_cols; ++i)
    {
        self->output_cols[i] = (char *) calloc(1, INITIAL_COL_SIZE *
                                               sizeof(char));
        // Make each col_ptrs pointer point to the beginning of the
        // column string
        self->col_ptrs[i] = self->output_cols[i];
        self->output_len[i] = INITIAL_COL_SIZE;
    }

    if (end == 0)
        RETURN(NO_ERROR); // don't read if end == 0

    self->state = START_LINE;

    // Loop until all of self->source has been read
    while (self->source_pos < self->source_len + 1)
    {
        if (self->source_pos == self->source_len || parse_newline)
            c = '\n';
        else
            c = self->source[self->source_pos];

        parse_newline = 0;

        switch (self->state)
        {
        case START_LINE:
            if (c == '\n')
                break;
            else if (c == '\r')
            {
                HANDLE_CR();
                break;
            }
            else if ((c == ' ' || c == '\t') && self->strip_whitespace_lines)
                break;
            else if (self->comment != 0 && c == self->comment)
            {
                // comment line; ignore
                self->state = COMMENT;
                break;
            }
            // initialize variables for the beginning of line parsing
            col = 0;
            BEGIN_FIELD();
            // parse in mode START_FIELD

        case START_FIELD:
            // strip whitespace before field begins
            if ((c == ' ' || c == '\t') && self->strip_whitespace_fields)
                break;
            else if (!self->strip_whitespace_lines && self->comment != 0 &&
                     c == self->comment)
            {
                // comment line, not caught earlier because of no stripping
                self->state = COMMENT;
                break;
            }
            else if (c == self->delimiter) // field ends before it begins
            {
                END_FIELD();
                BEGIN_FIELD();
                break;
            }
            else if (c == '\r')
            {
                HANDLE_CR();
                break;
            }
            else if (c == '\n')
            {
                if (self->strip_whitespace_lines)
                {
                    // Move on if the delimiter is whitespace, e.g.
                    // '1 2 3   '->['1','2','3']
                    if (self->delimiter == ' ' || self->delimiter == '\t')
                        ;
                    // Register an empty field if non-whitespace delimiter,
                    // e.g. '1,2, '->['1','2','']
                    else
                    {
                        END_FIELD();
                    }
                }

                else if (!self->strip_whitespace_lines)
                {
                    // In this case we don't want to left-strip the field,
                    // so we backtrack
                    int tmp = self->source_pos;
                    --self->source_pos;

                    while (self->source_pos >= 0 &&
                           self->source[self->source_pos] != self->delimiter
                           && self->source[self->source_pos] != '\n'
                           && self->source[self->source_pos] != '\r')
                    {
                        --self->source_pos;
                    }

                    // backtracked to line beginning
                    if (self->source_pos == -1 || self->source[self->source_pos] == '\n'
                        || self->source[self->source_pos] == '\r')
                    {
                        self->source_pos = tmp;
                    }
                    else
                    {
                        ++self->source_pos;

                        if (self->source_pos == tmp)
                            // no whitespace, just an empty field
                            ;

                        else
                            while (self->source_pos < tmp)
                            {
                                // append whitespace characters
                                PUSH(self->source[self->source_pos]);
                                ++self->source_pos;
                            }

                        if (col >= self->num_cols)
                            RETURN(TOO_MANY_COLS);
                        END_FIELD(); // whitespace counts as a field
                    }
                }

                END_LINE();
                self->state = START_LINE;
                break;
            }
            else if (c == self->quotechar) // start parsing quoted field
            {
                self->state = START_QUOTED_FIELD;
                break;
            }
            else
            {
                if (col >= self->num_cols)
                    RETURN(TOO_MANY_COLS);
                // Valid field character, parse again in FIELD mode
                self->state = FIELD;
            }

        case FIELD:
            if (self->comment != 0 && c == self->comment && whitespace && col == 0)
            {
                // No whitespace stripping, but the comment char is found
                // before any data, e.g. '  # a b c'
                self->state = COMMENT;
            }
            else if (c == self->delimiter)
            {
                // End of field, look for new field
                END_FIELD();
                BEGIN_FIELD();
            }
            else if (c == '\r')
            {
                HANDLE_CR();
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
                if (c != ' ' && c != '\t')
                    whitespace = 0; // field is not all whitespace
                PUSH(c);
            }
            break;

        case START_QUOTED_FIELD:
            if ((c == ' ' || c == '\t') && self->strip_whitespace_fields)
            {
                // ignore initial whitespace
                break;
            }
            else if (c == self->quotechar) // empty quotes
            {
                self->state = FIELD; // parse the rest of the field normally
            }
            else
            {
                // Valid field character, parse again in QUOTED_FIELD mode
                self->state = QUOTED_FIELD;
            }

        case QUOTED_FIELD_NEWLINE:
            if (self->state == QUOTED_FIELD)
                ; // fall through
            // Ignore initial whitespace if strip_whitespace_lines and
            // newlines regardless
            else if (((c == ' ' || c == '\t') && self->strip_whitespace_lines)
                     || c == '\n')
                break;
            else if (c == '\r')
            {
                HANDLE_CR();
                break;
            }
            else if (c == self->quotechar)
            {
                self->state = FIELD;
                break;
            }
            else
            {
                // Once data begins, parse it as a normal quoted field
                self->state = QUOTED_FIELD;
            }

        case QUOTED_FIELD:
            if (c == self->quotechar) // Parse rest of field normally, e.g. "ab"c
                self->state = FIELD;
            else if (c == '\n')
                self->state = QUOTED_FIELD_NEWLINE;
            else if (c == '\r')
            {
                HANDLE_CR();
            }
            else
            {
                PUSH(c);
            }
            break;

        case COMMENT:
            if (c == '\n')
            {
                self->state = START_LINE;
                if (!header)
                    end_comment(self);
            }
            else if (c == '\r')
            {
                HANDLE_CR();
            }
            else if (!header)
                push_comment(self, c);
            break; // keep looping until we find a newline

        case CARRIAGE_RETURN:
            self->state = old_state;
            --self->source_pos; // parse the newline in the old state
            if (c != '\n') // CR line terminator
            {
                --self->source_pos; // backtrack to the carriage return
                parse_newline = 1; // explicitly treat the CR as a newline
            }
            break;
        }

        ++self->source_pos;
    }

    RETURN(0);
}


// Lower-case a single C locale character
static inline int ascii_tolower(int c)
{
    if (c >= 'A' || c <= 'Z')
    {
        return c + ('a' - 'A');
    }

    return c;
}


static int ascii_strncasecmp(const char *str1, const char *str2, size_t n)
{
    int char1, char2;

    do
    {
        char1 = tolower(*(str1++));
        char2 = tolower(*(str2++));
        n--;
    } while (n && char1 != '\0' && char1 == char2);

    return (char1 - char2);
}


long str_to_long(tokenizer_t *self, char *str)
{
    char *tmp;
    long ret;
    errno = 0;
    ret = strtol(str, &tmp, 0);

    if (tmp == str || *tmp != '\0')
        self->code = CONVERSION_ERROR;
    else if (errno == ERANGE)
        self->code = OVERFLOW_ERROR;

    return ret;
}


double str_to_double(tokenizer_t *self, char *str)
{
    char *tmp;
    double val;
    errno = 0;

    if (self->use_fast_converter)
    {
        val = xstrtod(str, &tmp, '.', 'E', ',', 1);

        if (*tmp)
        {
            goto conversion_error;
        }
        else if (errno == ERANGE)
        {
            self->code = OVERFLOW_ERROR;
        }

        return val;
    }

    else
    {
        val = strtod(str, &tmp);

        if (errno == EINVAL || tmp == str || *tmp != '\0')
        {
            goto conversion_error;
        }
        else if (errno == ERANGE)
        {
            self->code = OVERFLOW_ERROR;
        }

        return val;
    }

conversion_error:
    // Handle inf and nan values for xstrtod and platforms whose strtod
    // doesn't support this
    val = 1.0;
    tmp = str;

    if (*tmp == '+')
    {
        tmp++;
    }
    else if (*tmp == '-')
    {
        tmp++;
        val = -1.0;
    }

    if (0 == ascii_strncasecmp(tmp, "nan", 3))
    {
        // Handle optional nan type specifier; this is ignored
        tmp += 3;
        val = NAN;
    }
    else if (0 == ascii_strncasecmp(tmp, "inf", 3))
    {
        tmp += 3;
        if (0 == ascii_strncasecmp(tmp, "inity", 5))
        {
            tmp += 5;
        }
        val *= INFINITY;
    }

    if (tmp == str || *tmp != '\0')
    {
        self->code = CONVERSION_ERROR;
        val = 0;
    }

    return val;
}

// ---------------------------------------------------------------------------
// Implementation of xstrtod

//
// strtod.c
//
// Convert string to double
//
// Copyright (C) 2002 Michael Ringgaard. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// -----------------------------------------------------------------------
// Modifications by Warren Weckesser, March 2011:
// * Rename strtod() to xstrtod().
// * Added decimal and sci arguments.
// * Skip trailing spaces.
// * Commented out the other functions.
// Modifications by Richard T Guy, August 2013:
// * Add tsep argument for thousands separator
// Modifications by Michael Mueller, August 2014:
// * Cache powers of 10 in memory to avoid rounding errors
// * Stop parsing decimals after 17 significant figures
//

double xstrtod(const char *str, char **endptr, char decimal,
               char sci, char tsep, int skip_trailing)
{
    double number;
    int exponent;
    int negative;
    char *p = (char *) str;
    int num_digits;
    int num_decimals;
    int max_digits = 17;
    int n;
    // Cache powers of 10 in memory
    static double e[] = {1., 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10,
                         1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20,
                         1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27, 1e28, 1e29, 1e30,
                         1e31, 1e32, 1e33, 1e34, 1e35, 1e36, 1e37, 1e38, 1e39, 1e40,
                         1e41, 1e42, 1e43, 1e44, 1e45, 1e46, 1e47, 1e48, 1e49, 1e50,
                         1e51, 1e52, 1e53, 1e54, 1e55, 1e56, 1e57, 1e58, 1e59, 1e60,
                         1e61, 1e62, 1e63, 1e64, 1e65, 1e66, 1e67, 1e68, 1e69, 1e70,
                         1e71, 1e72, 1e73, 1e74, 1e75, 1e76, 1e77, 1e78, 1e79, 1e80,
                         1e81, 1e82, 1e83, 1e84, 1e85, 1e86, 1e87, 1e88, 1e89, 1e90,
                         1e91, 1e92, 1e93, 1e94, 1e95, 1e96, 1e97, 1e98, 1e99, 1e100,
                         1e101, 1e102, 1e103, 1e104, 1e105, 1e106, 1e107, 1e108, 1e109, 1e110,
                         1e111, 1e112, 1e113, 1e114, 1e115, 1e116, 1e117, 1e118, 1e119, 1e120,
                         1e121, 1e122, 1e123, 1e124, 1e125, 1e126, 1e127, 1e128, 1e129, 1e130,
                         1e131, 1e132, 1e133, 1e134, 1e135, 1e136, 1e137, 1e138, 1e139, 1e140,
                         1e141, 1e142, 1e143, 1e144, 1e145, 1e146, 1e147, 1e148, 1e149, 1e150,
                         1e151, 1e152, 1e153, 1e154, 1e155, 1e156, 1e157, 1e158, 1e159, 1e160,
                         1e161, 1e162, 1e163, 1e164, 1e165, 1e166, 1e167, 1e168, 1e169, 1e170,
                         1e171, 1e172, 1e173, 1e174, 1e175, 1e176, 1e177, 1e178, 1e179, 1e180,
                         1e181, 1e182, 1e183, 1e184, 1e185, 1e186, 1e187, 1e188, 1e189, 1e190,
                         1e191, 1e192, 1e193, 1e194, 1e195, 1e196, 1e197, 1e198, 1e199, 1e200,
                         1e201, 1e202, 1e203, 1e204, 1e205, 1e206, 1e207, 1e208, 1e209, 1e210,
                         1e211, 1e212, 1e213, 1e214, 1e215, 1e216, 1e217, 1e218, 1e219, 1e220,
                         1e221, 1e222, 1e223, 1e224, 1e225, 1e226, 1e227, 1e228, 1e229, 1e230,
                         1e231, 1e232, 1e233, 1e234, 1e235, 1e236, 1e237, 1e238, 1e239, 1e240,
                         1e241, 1e242, 1e243, 1e244, 1e245, 1e246, 1e247, 1e248, 1e249, 1e250,
                         1e251, 1e252, 1e253, 1e254, 1e255, 1e256, 1e257, 1e258, 1e259, 1e260,
                         1e261, 1e262, 1e263, 1e264, 1e265, 1e266, 1e267, 1e268, 1e269, 1e270,
                         1e271, 1e272, 1e273, 1e274, 1e275, 1e276, 1e277, 1e278, 1e279, 1e280,
                         1e281, 1e282, 1e283, 1e284, 1e285, 1e286, 1e287, 1e288, 1e289, 1e290,
                         1e291, 1e292, 1e293, 1e294, 1e295, 1e296, 1e297, 1e298, 1e299, 1e300,
                         1e301, 1e302, 1e303, 1e304, 1e305, 1e306, 1e307, 1e308};
    errno = 0;

    // Skip leading whitespace
    while (isspace(*p)) p++;

    // Handle optional sign
    negative = 0;
    switch (*p)
    {
    case '-': negative = 1; // Fall through to increment position
    case '+': p++;
    }

    number = 0.;
    exponent = 0;
    num_digits = 0;
    num_decimals = 0;

    // Process string of digits
    while (isdigit(*p))
    {
        if (num_digits < max_digits)
        {
            number = number * 10. + (*p - '0');
            num_digits++;
        }
        else
            ++exponent;

        p++;
        p += (tsep != '\0' && *p == tsep);
    }

    // Process decimal part
    if (*p == decimal)
    {
        p++;

        while (num_digits < max_digits && isdigit(*p))
        {
            number = number * 10. + (*p - '0');
            p++;
            num_digits++;
            num_decimals++;
        }

        if (num_digits >= max_digits) // consume extra decimal digits
            while (isdigit(*p))
                ++p;

        exponent -= num_decimals;
    }

    if (num_digits == 0)
    {
        errno = ERANGE;
        return 0.0;
    }

    // Correct for sign
    if (negative) number = -number;

    // Process an exponent string
    if (toupper(*p) == toupper(sci))
    {
        // Handle optional sign
        negative = 0;
        switch (*++p)
        {
        case '-': negative = 1;   // Fall through to increment pos
        case '+': p++;
        }

        // Process string of digits
        n = 0;
        while (isdigit(*p))
        {
            n = n * 10 + (*p - '0');
            p++;
        }

        if (negative)
            exponent -= n;
        else
            exponent += n;
    }

    if (exponent > 308)
    {
        errno = ERANGE;
        return HUGE_VAL;
    }
    else if (exponent > 0)
        number *= e[exponent];
    else if (exponent < -308) // subnormal
    {
        if (exponent < -616) // prevent invalid array access
            number = 0.;
        number /= e[-308 - exponent];
        number /= e[308];
    }
    else
        number /= e[-exponent];

    if (number == HUGE_VAL || number == -HUGE_VAL)
        errno = ERANGE;

    if (skip_trailing) {
        // Skip trailing whitespace
        while (isspace(*p)) p++;
    }

    if (endptr) *endptr = p;
    return number;
}


void start_iteration(tokenizer_t *self, int col)
{
    // Begin looping over the column string with index col
    self->iter_col = col;
    // Start at the initial pointer position
    self->curr_pos = self->output_cols[col];
}


char *next_field(tokenizer_t *self, int *size)
{
    char *tmp = self->curr_pos;

    // pass through the entire field until reaching the delimiter
    while (*self->curr_pos != '\x00')
    ++self->curr_pos;

    ++self->curr_pos; // next field begins after the delimiter

    if (*tmp == '\x01') // empty field; this is a hack
    {
        if (size)
            *size = 0;
        return self->buf;
    }

    else
    {
        if (size)
            *size = self->curr_pos - tmp - 1;
        return tmp;
    }
}


char *get_line(char *ptr, int *len, int map_len)
{
    int pos = 0;

    while (pos < map_len)
    {
        if (ptr[pos] == '\r')
        {
            *len = pos;
            // Windows line break (\r\n)
            if (pos != map_len - 1 && ptr[pos + 1] == '\n')
                return ptr + pos + 2; // skip newline character
            else // Carriage return line break
                return ptr + pos + 1;
        }

        else if (ptr[pos] == '\n')
        {
            *len = pos;
            return ptr + pos + 1;
        }

        ++pos;
    }

    // done with input
    return 0;
}


void reset_comments(tokenizer_t *self)
{
    free(self->comment_lines);
    self->comment_pos = 0;
    self->comment_lines_len = INITIAL_COMMENT_LEN;
    self->comment_lines = (char *) malloc(INITIAL_COMMENT_LEN);
}
