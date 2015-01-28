import numpy as np

# Python 2 & 3 have different names for the __builtin__ module
try:
    import builtins 
except ImportError:
    import __builtin__ as builtins

def file_len(fname):
    """ Compute the number of all rows in fname

    Parameters 
    ----------
    fname : string 

    Returns 
    -------
    Nrows : int
 
    """
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    Nrows = i + 1
    return Nrows

def header_len(fname,header_char='#'):
    """ Compute the number of header rows in fname. 

    Parameters 
    ----------
    fname : string 

    header_char : str, optional
        string to be interpreted as a header line

    Returns 
    -------
    Nheader : int

    Notes 
    -----
    All empty lines that appear in header 
    will be included in the count. 

    """
    Nheader = 0
    with open(fname) as f:
        for i, l in enumerate(f):
            if ( (l[0:len(header_char)]==header_char) or (l=="\n") ):
                Nheader += 1
            else:
                break

    return Nheader

def get_header(fname, Nrows_header_total=None):
    """ Return the header as a list of strings, 
    one entry per header row. 

    Parameters 
    ----------
    fname : string 

    Nrows_header_total :  int, optional
        If the total number of header rows is not known in advance, 
        method will call `header_len` to determine Nrows_header_total. 

    Notes 
    -----
    Empty lines will be included in the returned header. 

    """

    if Nrows_header_total==None:
        Nrows_header_total = header_len(fname)

    output = []
    with open(fname) as f:
        for i in range(Nrows_header_total):
            line = f.readline().strip()
            output.append(line)

    return output


def read_catalog(fname, input_column_info, input_catalog_cuts=None,
    Nrows_total=None, Nrows_header_total=None):
    """

    Parameters 
    ----------
    fname : string 

    input_column_info : array_like 
        List whose entries are tuples providing information 
        about which columns of the ASCII data are desired. 

        There should be one tuple for each desired column. 
        Each tuple's first entry should give the column number 
        of the input ASCII file containing the desired data;
        the second entry should be any string 
        giving the desired name of the column 
        in the output data structure; 
        the third entry should be either 'float', 'int', or 'str', 
        depending on the data type of that column.

    input_catalog_cuts : array_like, optional 
        List whose entries are tuples providing information 
        about how to apply cuts to the columns. 
        There should be one tuple for each column 
        upon which a cut is being made. 
        Each tuple's first entry should give the column number 
        of the input ASCII file containing the desired data;
        the second entry should contain the desired lower bound 
        on the column (use None or float("-inf") if not applying a 
        lower bound cut); the third entry should contain 
        the desired upper bound on the column 
        (use None or float("inf") if not applying an upper bound cut). 

    Nrows_total : int, optional 
        If specified, the reader will only use 
        lines 0 through Nrows_total-1 of the ASCII file. 

    Nrows_header_total : int, optional 
        File lines 0 through Nrows_header_total-1 will be skipped. 
        Useful for non-trivial header formats that break the 
        simple `header_len` function. 

    Notes 
    -----
    Empty lines will automatically be skipped, including any number of  
    possibly empty lines at the end of the file. 

    """

    # If not provided, determine the number of header rows, and the total
    if Nrows_header_total==None:
        Nrows_header_total = header_len(fname)
    if Nrows_total==None:
        Nrows_total = file_len(fname)
    Ndata_rows = Nrows_total - Nrows_header_total

    # If applying cuts, make sure input_catalog_cuts contains a subset 
    # of columns in input_column_info. Append to input_column_info if not. 
    if input_catalog_cuts is not None:
        catalog_cuts = reformat_cuts(input_catalog_cuts) # Switches None cuts entries to float("inf")
        column_info = add_missing_columns(input_column_info, catalog_cuts)
    else:
        catalog_cuts=None
        column_info = input_column_info

    def get_builtin(name):
        """ Simple internal function to extract the appropriate method 
        to convert a string. 
        """
        t = getattr(builtins, name)
        if isinstance(t, type):
            return t
        raise ValueError(name)

    output = []

    with open(fname) as f:
        # Skip header
        for i in range(Nrows_header_total):
            f.readline()
        for i in range(Ndata_rows):
            parsed_line = f.readline().strip().split()
            if parsed_line != []:
                line = tuple(
                    get_builtin(column_info[idx][2])(parsed_line[column_info[idx][0]]) 
                    for idx in range(len(column_info)))
                if catalog_cuts == None:
                    output.append(line)
                else:
                    keep_row = catalogcut(parsed_line, column_info, catalog_cuts)
                    if keep_row == 1:
                        output.append(line)

    # Now bundle up the result into a structured numpy array
    def get_dt_string(type_list):
        """ Internal method to set the datatype of the 
        returned structured array
        """
        output = ''
        for entry in type_list:
            output = output+entry+'64,'
        return output[0:-1]

    colnames=np.array(column_info)[:,1].astype(str)
    dtt_string = get_dt_string(np.array(column_info)[:,2].astype(str))
    dt = np.dtype(dtt_string)
    output_array = np.array(output, dtype=dt)
    output_array.dtype.names = colnames

    return output_array


def single_variable_cut(x, xmin, xmax, 
    inclusion_convention=(False, False)):
    """ Method used to apply a cut to a single 
    variable in the catalog. 

    Parameters 
    ----------
    x : float or int 
        Value of the catalog upon which the cut is based 

    xmin : float or int 

    xmax : float or int 

    inclusion_convention : tuple of booleans 
        Used to determine whether upper and lower bounds 
        are inclusive or exclusive. False for exclusive, 
        True for inclusive; inclusion_convention[0] pertains 
        to lower bound, inclusion_convention[1] to upper bound. 

    Returns 
    -------
    output : boolean 
        Returns True if x passes the cut, False if not. 

    """

    if inclusion_convention==(True, True):
        return (x >= xmin) & (x <= xmax)
    elif inclusion_convention==(True, False):
        return (x >= xmin) & (x < xmax)
    elif inclusion_convention==(False, True):
        return (x > xmin) & (x <= xmax)
    else:
        return (x > xmin) & (x < xmax)

def catalogcut(data_row, column_info, catalogcuts):
    """ Determine whether the input catalog should be kept or discarded 
    based on a series of input cuts. 

    Parameters 
    ----------
    data_row : array_like 
        array of strings corresponding to a single row of the catalog table. 

    column_info : array_like 
        List of tuples provided as input to `read_catalog` to determine 
        which columns of the catalog are saved. 
        See docstring of `read_catalog` for details. 

    catalogcuts : array_like 
        List of tuples provided as input to `read_catalog` to determine 
        how cuts are applied to the catalog. 
        See docstring of `read_catalog` for details. 

    Returns 
    -------
    output : boolean 
        True if the row contains data for a catalog that passes the cuts, 
        False if not. 
    """
    output = True
    for cut in catalogcuts:
        # A row must pass ALL cuts in order to be included
        output *= single_variable_cut(float(data_row[cut[0]]), cut[1], cut[2])

    return output


def add_missing_columns(input_column_info, cuts):
    """ Augments input_column_info with columns used to 
    apply cuts to the catalog data, if applicable. 

    Parameters 
    ----------
    input_column_info : array_like 
        List of tuples provided as input to `read_catalog` to determine 
        which columns of the catalog are saved. 
        See docstring of `read_catalog` for details. 

    cuts : array_like 
        List of tuples provided as input to `read_catalog` to determine 
        how cuts are applied to the catalog. 
        See docstring of `read_catalog` for details. 

    Returns 
    -------
    output_column_info : array_like 
        List of tuples equal to input_column_info, plus possibly 
        tuple entries that appeared in cuts but not in input_column_info. 

    """
    output_column_info = input_column_info
    num_added_cols=0
    for cut in cuts:
        # the first entry of the cut tuple contains the 
        # relevant column number of the ASCII data
        colnum = cut[0]
        # Check to see whether the column was not included in 
        # the list of columns requested of the output data structure
        if colnum not in np.array(output_column_info)[:,0].astype(int):
            # Create a glaring, fake name for the column, 
            # assume it is a float, and add it to column_info
            num_added_cols += 1
            name = 'unlisted_col'+str(num_added_cols)+'_with_cut'
            datatype = 'float'
            new_tuple = (colnum, name, datatype)
            output_column_info.append(new_tuple)
        # Now the catalog reader will additionally include 
        # any columns used to make cuts that were omitted from 
        # the columns requested to be saved to an output data structure. 
        # This makes the coding simpler, and helps avoid the 
        # hidden-cut problem. 

    return output_column_info

def reformat_cuts(input_cuts):
    """ Convenience method that converts None entries of cuts list 
    to :math:`-\\infty` for lower bound, :math:`+\\infty` for upper bound. 
    """
    output_cuts = []
    for cut in input_cuts:
        cut = list(cut)
        if cut[1]==None:
            cut[1]=float("-inf")
        if cut[2]==None:
            cut[2]=float("inf")
        cut = tuple(cut)
        output_cuts.append(cut)
    return output_cuts









