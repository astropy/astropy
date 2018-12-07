# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from astropy.table.bst import BST


def get_tree(TreeType):
    b = TreeType([], [])
    for val in [5, 2, 9, 3, 4, 1, 6, 10, 8, 7]:
        b.add(val)
    return b


@pytest.fixture
def tree():
    return get_tree(BST)
    r'''
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


@pytest.fixture
def bst(tree):
    return tree


def test_bst_add(bst):
    root = bst.root
    assert root.data == [5]
    assert root.left.data == [2]
    assert root.right.data == [9]
    assert root.left.left.data == [1]
    assert root.left.right.data == [3]
    assert root.right.left.data == [6]
    assert root.right.right.data == [10]
    assert root.left.right.right.data == [4]
    assert root.right.left.right.data == [8]
    assert root.right.left.right.left.data == [7]


def test_bst_dimensions(bst):
    assert bst.size == 10
    assert bst.height == 4


def test_bst_find(tree):
    bst = tree
    for i in range(1, 11):
        node = bst.find(i)
        assert node == [i]
    assert bst.find(0) == []
    assert bst.find(11) == []
    assert bst.find('1') == []


def test_bst_traverse(bst):
    preord = [5, 2, 1, 3, 4, 9, 6, 8, 7, 10]
    inord = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    postord = [1, 4, 3, 2, 7, 8, 6, 10, 9, 5]
    traversals = {}
    for order in ('preorder', 'inorder', 'postorder'):
        traversals[order] = [x.key for x in bst.traverse(order)]
    assert traversals['preorder'] == preord
    assert traversals['inorder'] == inord
    assert traversals['postorder'] == postord


def test_bst_remove(bst):
    order = (6, 9, 1, 3, 7, 2, 10, 5, 4, 8)
    vals = set(range(1, 11))
    for i, val in enumerate(order):
        assert bst.remove(val) is True
        assert bst.is_valid()
        assert set([x.key for x in bst.traverse('inorder')]) == \
                vals.difference(order[:i+1])
        assert bst.size == 10 - i - 1
        assert bst.remove(-val) is False


def test_bst_duplicate(bst):
    bst.add(10, 11)
    assert bst.find(10) == [10, 11]
    assert bst.remove(10, data=10) is True
    assert bst.find(10) == [11]
    with pytest.raises(ValueError):
        bst.remove(10, data=30)  # invalid data
    assert bst.remove(10) is True
    assert bst.remove(10) is False


def test_bst_range(tree):
    bst = tree
    lst = bst.range_nodes(4, 8)
    assert sorted(x.key for x in lst) == [4, 5, 6, 7, 8]
    lst = bst.range_nodes(10, 11)
    assert [x.key for x in lst] == [10]
    lst = bst.range_nodes(11, 20)
    assert len(lst) == 0
