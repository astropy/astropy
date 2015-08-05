import operator
import numpy as np
from ..extern.six.moves import zip


class MaxValue(object):
    '''
    Represents an infinite value for purposes
    of tuple comparison.
    '''
    def __gt__(self, other):
        return True

    def __ge__(self, other):
        return True

    def __lt__(self, other):
        return False

    def __le__(self, other):
        return False

    def __repr__(self):
        return "MAX"

    __le__ = __lt__
    __ge__ = __gt__
    __str__ = __repr__


class MinValue(object):
    '''
    The opposite of MaxValue, i.e. a representation of
    negative infinity.
    '''
    def __lt__(self, other):
        return True

    def __le__(self, other):
        return True

    def __gt__(self, other):
        return False

    def __ge__(self, other):
        return False

    def __repr__(self):
        return "MIN"

    __le__ = __lt__
    __ge__ = __gt__
    __str__ = __repr__


class Epsilon(object):
    '''
    Represents the "next largest" version of a given value,
    so that for all valid comparisons we have
    x < y < Epsilon(y) < z whenever x < y < z and x, z are
    not Epsilon objects.

    Parameters
    ----------
    val : object
        Original value
    '''

    def __init__(self, val):
        self.val = val

    def __lt__(self, other):
        if self.val == other:
            return False
        return self.val < other

    def __gt__(self, other):
        if self.val == other:
            return True
        return self.val > other

    def __eq__(self, other):
        return False

    def __repr__(self):
        return repr(self.val) + " + epsilon"

class Node(object):
    '''
    An element in a binary search tree, containing
    a key, data, and references to children nodes and
    a parent node.

    Parameters
    ----------
    key : tuple
        Node key
    data : list or int
        Node data
    '''
    __lt__ = lambda x, y: x.key < y.key
    __le__ = lambda x, y: x.key <= y.key
    __eq__ = lambda x, y: x.key == y.key
    __ge__ = lambda x, y: x.key >= y.key
    __gt__ = lambda x, y: x.key > y.key
    __ne__ = lambda x, y: x.key != y.key

    # each node has a key and data list
    def __init__(self, key, data):
        self.key = key
        self.data = data if isinstance(data, list) else [data]
        self.left = None
        self.right = None
        self.parent = None

    def replace(self, child, new_child):
        '''
        Replace this node's child with a new child.
        '''
        if self.left is not None and self.left == child:
            self.set_left(new_child)
        elif self.right is not None and self.right == child:
            self.set_right(new_child)
        else:
            raise ValueError("Cannot call replace() on non-child")

    def remove(self, child):
        '''
        Remove the given child.
        '''
        self.replace(child, None)

    def set(self, other):
        '''
        Copy the given node.
        '''
        self.key = other.key
        self.data = other.data[:]

    def set_left(self, node):
        '''
        Set the left child to the given node.
        '''
        self.left = node
        if node is not None:
            node.parent = self

    def set_right(self, node):
        '''
        Set the right child to the given node.
        '''
        self.right = node
        if node is not None:
            node.parent = self

    def __str__(self):
        return str((self.key, self.data))

    def __repr__(self):
        return str(self)


class BST(object):
    '''
    A basic binary search tree in pure Python, used
    as an engine for indexing.

    Parameters
    ----------
    lines : `Table`
        Sorted columns of the original table as well as
        a final column consisting of the rows given by argsort()
    '''
    NodeClass = Node
    UNIQUE = False

    def __init__(self, lines=[]):
        self.root = None
        self.size = 0
        for row in lines:
            row = tuple(row)
            key, data = row[:-1], row[-1]
            self.add(key, data)

    def add(self, key, data=None):
        '''
        Add a key, data pair.
        '''
        if data is None:
            data = key

        self.size += 1
        node = self.NodeClass(key, data)
        curr_node = self.root
        if curr_node is None:
            self.root = node
            return
        while True:
            if node < curr_node:
                if curr_node.left is None:
                    curr_node.set_left(node)
                    break
                curr_node = curr_node.left
            elif node > curr_node:
                if curr_node.right is None:
                    curr_node.set_right(node)
                    break
                curr_node = curr_node.right
            elif self.UNIQUE:
                raise ValueError("Cannot insert non-unique value")
            else: # add data to node
                curr_node.data.extend(node.data)
                curr_node.data = sorted(curr_node.data)
                return

    def find(self, key):
        '''
        Return all data values corresponding to a given key.

        Parameters
        ----------
        key : tuple
            Input key

        Returns
        -------
        data_vals : list
            List of rows corresponding to the input key
        '''
        node = self.find_node(key)
        return node.data if node is not None else []

    def find_node(self, key):
        '''
        Find the node associated with the given key.
        '''
        if self.root is None:
            return None
        return self._find_recursive(key, self.root)

    def shift_left(self, row):
        '''
        Decrement all rows larger than the given row.
        '''
        for node in self.traverse():
            node.data = [x - 1 if x > row else x for x in node.data]

    def shift_right(self, row):
        '''
        Increment all rows greater than or equal to the given row.
        '''
        for node in self.traverse():
            node.data = [x + 1 if x >= row else x for x in node.data]

    def _find_recursive(self, key, node):
        try:
            if key == node.key:
                return node
            elif key > node.key:
                if node.right is None:
                    return None
                return self._find_recursive(key, node.right)
            else:
                if node.left is None:
                    return None
                return self._find_recursive(key, node.left)
        except TypeError: # wrong key type
            return None

    def traverse(self, order='inorder'):
        '''
        Return nodes of the BST in the given order.

        Parameters
        ----------
        order : str
            The order in which to recursively search the BST.
            Possible values are:
            "preorder": current node, left subtree, right subtree
            "inorder": left subtree, current node, right subtree
            "postorder": left subtree, right subtree, current node
        '''
        if order == 'preorder':
            return self._preorder(self.root, [])
        elif order == 'inorder':
            return self._inorder(self.root, [])
        elif order == 'postorder':
            return self._postorder(self.root, [])
        raise ValueError("Invalid traversal method: \"{0}\"".format(order))

    def items(self):
        '''
        Return BST items in order as (key, data) pairs.
        '''
        return [(x.key, x.data) for x in self.traverse()]

    def sort(self):
        '''
        Return BST rows sorted by key values.
        '''
        return [x for node in self.traverse() for x in node.data]

    def _preorder(self, node, lst):
        if node is None:
            return lst
        lst.append(node)
        self._preorder(node.left, lst)
        self._preorder(node.right, lst)
        return lst

    def _inorder(self, node, lst):
        if node is None:
            return lst
        self._inorder(node.left, lst)
        lst.append(node)
        self._inorder(node.right, lst)
        return lst

    def _postorder(self, node, lst):
        if node is None:
            return lst
        self._postorder(node.left, lst)
        self._postorder(node.right, lst)
        lst.append(node)
        return lst

    def _substitute(self, node, new_node):
        if node is self.root:
            self.root = new_node
        else:
            node.parent.replace(node, new_node)

    def remove(self, key, data=None):
        '''
        Remove data corresponding to the given key.

        Parameters
        ----------
        key : tuple
            The key to remove
        data : int or None
            If None, remove the node corresponding to the given key.
            If not None, remove only the given data value from the node.

        Returns
        -------
        successful : bool
            True if removal was successful, false otherwise
        '''
        node = self.find_node(key)
        if node is None:
            return False
        if data is not None:
            if data not in node.data:
                raise ValueError("Data does not belong to correct node")
            elif len(node.data) > 1:
                node.data.remove(data)
                return True
        if node.left is None and node.right is None:
            self._substitute(node, None)
        elif node.left is None and node.right is not None:
            self._substitute(node, node.right)
        elif node.right is None and node.left is not None:
            self._substitute(node, node.left)
        else:
            # find largest element of left subtree
            curr_node = node.left
            while curr_node.right is not None:
                curr_node = curr_node.right
            self._substitute(curr_node, curr_node.left)
            node.set(curr_node)
        self.size -= 1
        return True

    def is_valid(self):
        '''
        Returns whether this is a valid BST.
        '''
        return self._is_valid(self.root)

    def _is_valid(self, node):
        if node is None:
            return True
        return (node.left is None or node.left <= node) and \
            (node.right is None or node.right >= node) and \
            self._is_valid(node.left) and self._is_valid(node.right)

    def range(self, lower, upper, bounds=(True, True)):
        '''
        Return all nodes with keys in the given range.

        Parameters
        ----------
        lower : tuple
            Lower bound
        upper : tuple
            Upper bound
        bounds : tuple (x, y) of bools
            Indicates whether the search should be inclusive or
            exclusive with respect to the endpoints. The first
            argument x corresponds to an inclusive lower bound,
            and the second argument y to an inclusive upper bound.
        '''
        nodes = self.range_nodes(lower, upper, bounds)
        return [x for node in nodes for x in node.data]

    def range_nodes(self, lower, upper, bounds=(True, True)):
        '''
        Return nodes in the given range.
        '''
        if self.root is None:
            return []
        # ops are <= or <
        ops = tuple([operator.le if x else operator.lt for x in bounds])
        return self._range(lower, upper, ops, self.root, [])

    def same_prefix(self, val):
        '''
        Assuming the given value has smaller length than keys, return
        nodes whose keys have this value as a prefix.
        '''
        if self.root is None:
            return []
        nodes = self._same_prefix(val, self.root, [])
        return [x for node in nodes for x in node.data]

    def _range(self, lower, upper, ops, node, lst):
        op1, op2 = ops
        if op1(lower, node.key) and op2(node.key, upper):
            lst.append(node)
        if node.key < upper and node.right is not None:
            self._range(lower, upper, ops, node.right, lst)
        if node.key > lower and node.left is not None:
            self._range(lower, upper, ops, node.left, lst)
        return lst

    def _same_prefix(self, val, node, lst):
        prefix = node.key[:len(val)]
        if prefix == val:
            lst.append(node)
        if prefix <= val and node.right is not None:
            self._same_prefix(val, node.right, lst)
        if prefix >= val and node.left is not None:
            self._same_prefix(val, node.left, lst)
        return lst

    def __str__(self):
        if self.root is None:
            return 'Empty'
        return self._print(self.root, 0)

    def __repr__(self):
        return str(self)

    def _print(self, node, level):
        line = '\t'*level + str(node) + '\n'
        if node.left is not None:
            line += self._print(node.left, level + 1)
        if node.right is not None:
            line += self._print(node.right, level + 1)
        return line

    def height(self):
        '''
        Return the BST height.
        '''
        return self._height(self.root)

    def _height(self, node):
        if node is None:
            return -1
        return max(self._height(node.left),
                   self._height(node.right)) + 1

    def replace_rows(self, row_map):
        '''
        Replace all rows with the values they map to in the
        given dictionary. Any rows not present as keys in
        the dictionary will have their nodes deleted.

        Parameters
        ----------
        row_map : dict
            Mapping of row numbers to new row numbers
        '''
        for key, data in self.items():
            data[:] = [row_map[x] for x in data if x in row_map]


class FastBase(object):
    '''
    A fast binary search tree implementation for indexing,
    using the bintrees library.

    Parameters
    ----------
    lines : `Table`
        Sorted columns of the original table as well as
        a final column consisting of the rows given by argsort()
    '''
    def __init__(self, lines):
        self.data = self.engine()
        for row in np.array(lines):
            row = tuple(row)
            self.add(row[:-1], row[-1])

    def add(self, key, val):
        '''
        Add a key, value pair.
        '''
        rows = self.data.set_default(key, [])
        rows.insert(np.searchsorted(rows, val), val)

    def find(self, key):
        '''
        Find rows corresponding to the given key.
        '''
        return self.data.get(key, [])

    def remove(self, key, data=None):
        '''
        Remove data from the given key.
        '''
        node = self.data.get(key, None)
        if node is None or len(node) == 0:
            return False
        if data is None:
            self.data.pop(key)
            return True
        if data not in node:
            if len(node) == 0:
                return False
            raise ValueError("Data does not belong to correct node")
        node.remove(data)
        return True


    def shift_left(self, row):
        '''
        Decrement rows larger than the given row.
        '''
        for key, node in self.data.items():
            self.data[key] = [x - 1 if x > row else x for x in node]

    def shift_right(self, row):
        '''
        Increment rows greater than or equal to the given row.
        '''
        for key, node in self.data.items():
            self.data[key] = [x + 1 if x >= row else x for x in node]

    def traverse(self):
        '''
        Return all nodes in this BST.
        '''
        l = []
        for key, data in self.data.items():
            n = Node(key, key)
            n.data = data
            l.append(n)
        return l

    def items(self):
        '''
        Return a list of key, data tuples.
        '''
        return [x for x in self.data.items() if len(x[1]) > 0]

    def sort(self):
        '''
        Return a list of rows in order sorted by key.
        '''
        return [x for node in self.data.values() for x in node]

    def range(self, lower, upper, bounds=(True, True)):
        '''
        Return row values in the given range.
        '''
        # we need Epsilon since bintrees searches for
        # lower <= key < upper, while we might want lower <= key <= upper
        # or similar
        if not bounds[0]: # lower < key
            lower = Epsilon(lower)
        if bounds[1]: # key <= upper
            upper = Epsilon(upper)
        l = [v for v in self.data.value_slice(lower, upper)]
        return [x for sublist in l for x in sublist]

    def replace_rows(self, row_map):
        '''
        Replace rows with the values in row_map.
        '''
        for data in self.data.values():
            data[:] = [row_map[x] for x in data if x in row_map]

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return str(self)

try:
    # bintrees is an optional dependency
    from bintrees import FastBinaryTree, FastRBTree

    class FastBST(FastBase):
        engine = FastBinaryTree

    class FastRBT(FastBase):
        engine = FastRBTree

except ImportError:
    FastBST = BST
    FastRBT = BST
