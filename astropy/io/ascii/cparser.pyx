# Licensed under a 3-clause BSD style license - see LICENSE.rst

import six
import numpy as np
cimport numpy as np
from numpy import ma
from ...utils.data import get_readable_fileobj

cdef extern from "src/tokenizer.h":
	ctypedef enum tokenizer_state:
		START_LINE
		START_FIELD
		START_QUOTED_FIELD
		FIELD
		QUOTED_FIELD
		QUOTED_FIELD_NEWLINE
		COMMENT

	ctypedef enum err_code:
		NO_ERROR
		INVALID_LINE
		TOO_MANY_COLS

	ctypedef struct tokenizer_t:
		char *source		# single string containing all of the input
		int source_len 		# length of the input
		int source_pos		# current index in source for tokenization
		char delimiter		# delimiter character
		char comment		# comment character
		char quotechar		# quote character
		char *header_output # string containing header data
		char **output_cols	# array of output strings for each column
		char **col_ptrs     # array of pointers to current output position for each col
		int *output_len		# length of each output column string
		int header_len      # length of the header output string
		int num_cols		# number of table columns
		int num_rows		# number of table rows
		tokenizer_state state   # current state of the tokenizer
		err_code code		# represents the latest error that has occurred
		# Example input/output
		# --------------------
		# source: "A,B,C\n10,5.,6\n1,2,3"
		# output_cols: ["A101", "B5.2", "C6 3"]
		# row_positions: [0, 1, 3]

	tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar)
	void delete_tokenizer(tokenizer_t *tokenizer)
	int tokenize(tokenizer_t *self, int start, int end, int header, int *use_cols)
	int int_size()

class CParserError(Exception):
	"""
	An instance of this class is thrown when an error occurs
	during C parsing.
	"""

ERR_CODES = dict(enumerate([
	"no error",
	"invalid line supplied",
	lambda line: "too many columns found in line {} of data".format(line)
	]))

cdef class CParser:
	"""
	A fast Cython parser class which uses underlying C code
	for tokenization.
	"""

	cdef:
		tokenizer_t *tokenizer
		object source
		object header_start
		int data_start
		int data_end
		object data_end_obj
		object include_names
		object exclude_names
		object fill_values
		object fill_include_names
		object fill_exclude_names
		object fill_names
		np.ndarray use_cols

	cdef public:
		int width
		object names

	def __cinit__(self, source,
				  delimiter=',',
				  comment=None,
				  quotechar='"',
				  header_start=0,
				  data_start=1,
				  data_end=None,
				  names=None,
				  include_names=None,
				  exclude_names=None,
				  fill_values=('', '0'),
				  fill_include_names=None,
				  fill_exclude_names=None):

		self.tokenizer = create_tokenizer(ord(delimiter), ord(comment), ord(quotechar))
		self.source = None
		self.setup_tokenizer(source)
		self.header_start = header_start
		self.data_start = data_start
		self.data_end = -1 # keep reading data until the end
		if data_end is not None and data_end >= 0:
			self.data_end = data_end
		self.data_end_obj = data_end
		self.names = names
		self.include_names = include_names
		self.exclude_names = exclude_names
		if len(fill_values) > 0 and isinstance(fill_values[0], six.string_types):
			self.fill_values = [fill_values]
		else:
			self.fill_values = fill_values
		try:
			# Create a dict with the values to be replaced as keys
			self.fill_values = dict([(l[0], l[1:]) for l in self.fill_values])
		except IndexError:
			raise ValueError("Format of fill_values must be "
							 "(<bad>, <fill>, <optional col1>, ...)")
		self.fill_include_names = fill_include_names
		self.fill_exclude_names = fill_exclude_names
	
	def __dealloc__(self):
		if self.tokenizer:
			delete_tokenizer(self.tokenizer)

	cdef raise_error(self, msg):
		err_msg = ERR_CODES.get(self.tokenizer.code, "unknown error")
		if callable(err_msg):
			err_msg = err_msg(self.tokenizer.num_rows + 1)
		raise CParserError("{}: {}".format(msg, err_msg))

	cdef setup_tokenizer(self, source):
		cdef char *src

		if isinstance(source, six.string_types) or hasattr(source, 'read'):
			#todo: handle the case where source is the actual data (includes newline)
			with get_readable_fileobj(source) as file_obj:
				source = file_obj.read()
		else:
			try:
				source = '\n'.join(source) # iterable sequence of lines
			except TypeError:
				raise TypeError('Input "table" must be a file-like object, a string (filename'
			  			   'or data), or an iterable')
		# Create a reference to the Python object so its char * pointer remains valid
		self.source = source + '\n' # add newline to simplify handling last line of data
		src = self.source
		self.tokenizer.source = src
		self.tokenizer.source_len = len(self.source)

	def read_header(self):
		if self.names:
			self.width = len(self.names)
		# header_start is a valid line number
		elif self.header_start is not None and self.header_start >= 0:
			if tokenize(self.tokenizer, self.header_start, -1, 1, <int *> 0) != 0:
				self.raise_error("an error occurred while tokenizing the header line")
			self.names = []
			name = ''
			for i in range(self.tokenizer.header_len):
				c = self.tokenizer.header_output[i]
				if not c:
					if name:
						self.names.append(name.replace('\x01', ''))
						name = ''
					else:
						break # end of string
				else:
					name += chr(c)
			self.width = len(self.names)
		else:
			# Get number of columns from first data row
			if tokenize(self.tokenizer, 0, -1, 1, <int *> 0) != 0:
				self.raise_error("an error occurred while tokenizing the first line of data")
			self.width = 0
			for i in range(self.tokenizer.header_len):
				if not self.tokenizer.header_output[i]:
					if i > 0 and self.tokenizer.header_output[i - 1]:
						self.width += 1
					else:
						break
			self.names = ['col{}'.format(i + 1) for i in range(self.width)] # auto-generate names

		size = int_size()
		dtype = np.int16 #TODO: maybe find a better way to do this?
		if size == 64:
			dtype = np.int64
		elif size == 32:
			dtype = np.int32
		self.use_cols = np.ones(self.width, dtype)
		if self.include_names is not None:
			for i, name in enumerate(self.names):
				if name not in self.include_names:
					self.use_cols[i] = 0
		if self.exclude_names is not None:
			for name in self.exclude_names:
				try:
					self.use_cols[self.names.index(name)] = 0
				except ValueError: # supplied name is invalid, ignore
					continue

		self.names = [self.names[i] for i, should_use in enumerate(self.use_cols) if should_use]
		self.width = len(self.names)
		self.tokenizer.num_cols = self.width
			
	def read(self):
		if tokenize(self.tokenizer, self.data_start, self.data_end, 0,
					<int *> self.use_cols.data) != 0:
			self.raise_error("an error occurred while tokenizing data")
		else:
			self._set_fill_names()
			return self._convert_data()

	cdef _set_fill_names(self):
		self.fill_names = set(self.names)
		if self.fill_include_names is not None:
			self.fill_names.intersection_update(self.fill_include_names)
		if self.fill_exclude_names is not None:
			self.fill_names.difference_update(self.fill_exclude_names)

	cdef _convert_data(self):
		cols = {}
		cdef int row
		cdef int num_rows = self.tokenizer.num_rows
		if self.data_end_obj is not None and self.data_end_obj < 0:
			num_rows += self.data_end_obj

		for i in range(self.tokenizer.num_cols):
			cols[self.names[i]] = np.empty(num_rows, dtype='|S{}'.format(self.tokenizer.output_len[i]))
			el = ''
			row = 0
			masked = False

			for j in range(self.tokenizer.output_len[i]):
				if row >= num_rows:
					break
				c = self.tokenizer.output_cols[i][j]
				if not c:
					if not el:
						break
					if el == '\x01':
						el = ''
					if el in self.fill_values: #TODO: tighten this bit up clarity-wise
						new_val = str(self.fill_values[el][0])
						if (len(self.fill_values[el]) > 1 and self.names[i] in self.fill_values[el][1:]) or \
						   (len(self.fill_values[el]) == 1 and self.names[i] in self.fill_names):
							if not masked:
								masked = True
								# change from ndarray to MaskedArray
								cols[self.names[i]] = cols[self.names[i]].view(ma.MaskedArray)
								# by default, values are not masked
								cols[self.names[i]].mask = np.zeros(num_rows)
							cols[self.names[i]][row] = new_val
							cols[self.names[i]].mask[row] = True
						else:
							cols[self.names[i]][row] = el
					else:
						cols[self.names[i]][row] = el
					el = ''
					row += 1
				else:
					el += chr(c)

		for name in self.names:
			for dtype in (np.int_, np.float_):
				try:
					cols[name] = cols[name].astype(dtype)
					break
				except ValueError:
					continue
		return cols
