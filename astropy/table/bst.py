class Node:
    __lt__ = lambda x, y: x.key < y.key
    __le__ = lambda x, y: x.key <= y.key
    __eq__ = lambda x, y: x.key == y.key
    __ge__ = lambda x, y: x.key >= y.key
    __gt__ = lambda x, y: x.key > y.key
    __ne__ = lambda x, y: x.key != y.key

    # each node has a key and data list
    def __init__(self, key, data=None,
                 left=None, right=None, parent=None):
        self.key = key
        self.data = [data] if data is not None else [key]
        self.left = left
        self.right = right
        self.parent = parent

    def replace(self, child, new_child):
        if self.left is not None and self.left == child:
            self.left = new_child
            if new_child is not None:
                new_child.parent = self
        elif self.right is not None and self.right == child:
            self.right = new_child
            if new_child is not None:
                new_child.parent = self
        else:
            raise ValueError("Cannot call replace() on non-child")

    def remove(self, child):
        self.replace(child, None)

    def set(self, other):
        self.key = other.key
        self.data = other.data[:]

class BST:
    NodeClass = Node
    UNIQUE = False

    def __init__(self, data=[]):
        self.root = None
        self.size = 0
        for d in data:
            self.add(*d)

    def add(self, *data):
        self.size += 1
        node = self.NodeClass(*data)
        curr_node = self.root
        if curr_node is None:
            self.root = node
            return
        while True:
            if node < curr_node:
                if curr_node.left is None:
                    curr_node.left = node
                    node.parent = curr_node
                    break
                curr_node = curr_node.left
            elif node > curr_node:
                if curr_node.right is None:
                    curr_node.right = node
                    node.parent = curr_node
                    break
                curr_node = curr_node.right
            elif self.UNIQUE:
                raise ValueError("Cannot insert non-unique value")
            else: # add data to node
                curr_node.data.extend(node.data)
                curr_node.data = sorted(curr_node.data) ##TODO: speed up
                return
        self.balance()

    def balance(self):
        pass

    def find(self, key):
        if self.root is None:
            return None
        return self._find_recursive(key, self.root)

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

    def traverse(self, order):
        if order == 'preorder':
            return self._preorder(self.root, [])
        elif order == 'inorder':
            return self._inorder(self.root, [])
        elif order == 'postorder':
            return self._postorder(self.root, [])
        raise ValueError("Invalid traversal method: \"{0}\"".format(order))
        ##TODO: find out why inorder is so slow

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
        node = self.find(key)
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
        if self.root is None:
            return []
        return self._range(lower, upper, self.root, [])

    def _range(self, lower, upper, node, lst):
        if lower <= node.key <= upper:
            lst.append(node)
        if node.key < upper and node.right is not None:
            self._range(lower, upper, node.right, lst)
        if node.key > lower and node.left is not None:
            self._range(lower, upper, node.left, lst)
        return lst
