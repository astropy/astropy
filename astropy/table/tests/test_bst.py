from ..bst import BST

bst = BST()
for val in [5, 2, 9, 3, 4, 1, 6, 10, 8, 7]:
        bst.add(val)
'''
         5
       /   \
      2     9
     / \   / \
    1   3 6  10
         \ \
         4  8
           /
          7
'''

def test_bst_add():
    root = bst.root
    assert root.data == 5
    assert root.left.data == 2
    assert root.right.data == 9
    assert root.left.left.data == 1
    assert root.left.right.data == 3
    assert root.right.left.data == 6
    assert root.right.right.data == 10
    assert root.left.right.right.data == 4
    assert root.right.left.right.data == 8
    assert root.right.left.right.left.data == 7

def test_bst_size():
    assert bst.size == 10

def test_bst_find():
    for i in range(1, 11):
        node = bst.find(i)
        assert node is not None
        assert node == i
    assert bst.find(0) is None
    assert bst.find(11) is None
    assert bst.find('1') is None

def test_bst_traverse():
    preord = [5, 2, 1, 3, 4, 9, 6, 8, 7, 10]
    inord = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    postord = [1, 4, 3, 2, 7, 8, 6, 10, 9, 5]
    traversals = {}
    for order in ('preorder', 'inorder', 'postorder'):
        traversals[order] = [x.key for x in bst.traverse(order)]
    assert traversals['preorder'] == preord
    assert traversals['inorder'] == inord
    assert traversals['postorder'] == postord
