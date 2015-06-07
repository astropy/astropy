class Node:
    def __init__(self, key, data=None, left=None, right=None):
        self.key = key
        self.data = data if data is not None else key
        self.left = left
        self.right = right

    def __cmp__(self, other):
        if isinstance(other, Node):
            return self.key.__cmp__(other.key)
        return self.key.__cmp__(other)

class BST:
    NodeClass = Node

    def __init__(self):
        self.root = None
        self.size = 0

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
                    break
                curr_node = curr_node.left
            elif node > curr_node:
                if curr_node.right is None:
                    curr_node.right = node
                    break
                curr_node = curr_node.right
        self.balance()

    def balance(self):
        pass

    def find(self, key):
        return self._find_recursive(key, self.root)

    def _find_recursive(self, key, node):
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

    def traverse(self, order):
        if order == 'preorder':
            return self._preorder(self.root, [])
        elif order == 'inorder':
            return self._inorder(self.root, [])
        elif order == 'postorder':
            return self._postorder(self.root, [])
        raise ValueError("Invalid traversal method: \"{0}\"".format(order))

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
