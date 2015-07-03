BLACK = 0
RED = 1

class Epsilon(object):
    # Represents the "next largest" version of a given value,
    # so for all valid comparisons we have
    # x < y < Epsilon(y) < z wherever x < y < z and x, z are
    # not Epsilon objects

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

class Node(object):
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
        if self.left is not None and self.left == child:
            self.set_left(new_child)
        elif self.right is not None and self.right == child:
            self.set_right(new_child)
        else:
            raise ValueError("Cannot call replace() on non-child")

    def remove(self, child):
        self.replace(child, None)

    def set(self, other):
        self.key = other.key
        self.data = other.data[:]

    def set_left(self, node):
        self.left = node
        if node is not None:
            node.parent = self

    def set_right(self, node):
        self.right = node
        if node is not None:
            node.parent = self

    def __str__(self):
        return str((self.key, self.data))

    def __repr__(self):
        return str(self)

class RedBlackNode(Node):
    def __init__(self, key, data):
        super(RedBlackNode, self).__init__(key, data)
        self.color = RED

    def __str__(self):
        return str((self.key, self.data, 'BLACK' if self.color == BLACK else 'RED'))
    def __repr__(self):
        return str(self)

class BST(object):
    NodeClass = Node
    UNIQUE = False

    def __init__(self, lines={}):
        self.root = None
        self.size = 0
        for key, data in lines.items():
            self.add(key, data)

    def add(self, key, data=None):
        if data is None:
            data = key

        self.size += 1
        node = self.NodeClass(key, data)
        curr_node = self.root
        if curr_node is None:
            self.root = node
            node.color = BLACK
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
                curr_node.data = sorted(curr_node.data) ##TODO: speed up
                return
        self.balance(node)

    def balance(self, node):
        pass

    def find(self, key):
        node = self.find_node(key)
        return node.data if node is not None else []

    def find_node(self, key):
        if self.root is None:
            return None
        return self._find_recursive(key, self.root)

    def reorder(self, row):
        for node in self.traverse():
            node.data = [x - 1 if x > row else x for x in node.data]

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
        if order == 'preorder':
            return self._preorder(self.root, [])
        elif order == 'inorder':
            return self._inorder(self.root, [])
        elif order == 'postorder':
            return self._postorder(self.root, [])
        raise ValueError("Invalid traversal method: \"{0}\"".format(order))

    def sort(self):
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
        # returns True if successfully removed, False otherwise
        # if data is not None, remove the entire node if only this data
        # is present, otherwise just pop off data from the node
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
        return self._is_valid(self.root)

    def _is_valid(self, node):
        if node is None:
            return True
        return (node.left is None or node.left <= node) and \
            (node.right is None or node.right >= node) and \
            self._is_valid(node.left) and self._is_valid(node.right)

    def range(self, lower, upper):
        # return all nodes with keys in (inclusive) range [lower, upper]
        nodes = self.range_nodes(lower, upper)
        return [x for node in nodes for x in node.data]

    def range_nodes(self, lower, upper):
        if self.root is None:
            return []
        return self._range(lower, upper, self.root, [])

    def same_prefix(self, val):
        # assuming val has smaller length than keys, return
        # nodes whose keys have val as a prefix
        if self.root is None:
            return []
        nodes = self._same_prefix(val, self.root, [])
        return [x for node in nodes for x in node.data]

    def _range(self, lower, upper, node, lst):
        if lower <= node.key <= upper:
            lst.append(node)
        if node.key < upper and node.right is not None:
            self._range(lower, upper, node.right, lst)
        if node.key > lower and node.left is not None:
            self._range(lower, upper, node.left, lst)
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

    def nodes(self):
        # for debugging
        return [(x.key, x.data) for x in self.traverse('inorder')]

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
        return self._height(self.root)

    def _height(self, node):
        if node is None:
            return -1
        return max(self._height(node.left),
                   self._height(node.right)) + 1

    def replace_rows(self, row_map):
        for key, data in self.nodes():
            data[:] = [row_map[x] for x in data if x in row_map]

    def rotate_left(self, node):
        parent = node.parent
        subtree = node.right.left
        new_node = node.right
        node.set_right(subtree)
        new_node.set_left(node)

        if parent is not None:
            parent.replace(node, new_node)
        else:
            self.root = new_node
            new_node.parent = None

    def rotate_right(self, node):
        parent = node.parent
        subtree = node.left.right
        new_node = node.left
        node.set_left(subtree)
        new_node.set_right(node)

        if parent is not None:
            parent.replace(node, new_node)
        else:
            self.root = new_node
            new_node.parent = None

class RedBlackTree(BST):
    NodeClass = RedBlackNode

    def balance(self, node):
        node.color = RED
        parent = node.parent
        if parent is None:
            node.color = BLACK
            if node is not self.root:
                raise ValueError("Non-root node has no parent")
            return
        gp = parent.parent
        uncle = None
        if gp is not None:
            uncle = gp.right if gp.left is parent else gp.left
        elif parent.color == BLACK:
            pass
        elif uncle is not None and uncle.color == RED:
            # parent and uncle are red
            parent.color = BLACK
            uncle.color = BLACK
            gp.color == RED
            self.balance(gp)
        else: # parent is red, uncle is black (or None)
            # note: gp is not None because the root is black, handled above
            if node is parent.right and parent is gp.left:
                self.rotate_left(parent)
                node = node.left
            elif node is parent.left and parent is gp.right:
                self.rotate_right(parent)
                node = node.right
            # left/left or right/right
            node.parent.color = BLACK
            node.parent.parent.color = RED
            if node is node.parent.left:
                self.rotate_right(node.parent.parent)
            else:
                self.rotate_left(node.parent.parent)

    def is_valid(self):
        if self.root is None:
            return True
        return super(RedBlackTree, self).is_valid() and \
            self.root.color == BLACK and self._rbt_check(self.root)

    def _is_black(self, node):
        return node is None or node.color == BLACK

    def _rbt_check(self, node):
        ##TODO: check black-height property
        if node.color == RED:
            if not (self._is_black(node.left) and self._is_black(node.right)):
                return False
        if node.left is not None:
            return self._rbt_check(node.left)
        if node.right is not None:
            return self._rbt_check(node.right)
        return True

    ##TODO: write remove method

class FastBase(object):
    def __init__(self, lines):
        self.data = self.engine(lines)

    def add(self, key, val):
        self.data.set_default(key, []).append(val)

    def find(self, key):
        return self.data.get(key, None)

    def remove(self, key, data=None):
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

    def reorder(self, row):
        for key, node in self.data.items():
            self.data[key] = [x - 1 if x > row else x for x in node]

    def traverse(self):
        l = []
        for key, data in self.data.items():
            if len(data) > 0:
                n = Node(key, key)
                n.data = data
                l.append(n)
        return l

    def sort(self):
        return [x for node in self.traverse() for x in node.data]

    def range(self, lower, upper):
        # we need Epsilon(upper) since bintrees searches for
        # lower <= key < upper, while we want lower <= key <= upper
        l = [v for v in self.data.value_slice(lower, Epsilon(upper))]
        return [x for sublist in l for x in sublist]

    def replace_rows(self, row_map):
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
    FastRBT = BST ##TODO: change to RedBlackTree when completed
