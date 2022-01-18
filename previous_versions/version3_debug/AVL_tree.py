class AVLTreeNode:
    def __init__(self, key, val=None, bf=0):
        self.key = key
        self.val = val
        self.left = None  # Left child pointer
        self.right = None  # Right child pointer
        self.parent = None  # + Parent of this node for easier rotations
        self.bf = bf  # + Balance factor of current node
        self.height = 1  # + Height of subtree rooted at this node important to update bf's after re-balancing

    def __str__(self):
        return f'({self.key}, {self.val}, {self.height})'


class AVLTree:
    def __init__(self):
        self.root = None  # root of tree
        self.count = 0

    def _get_height(self, root: AVLTreeNode):
        if not root:  # empty tree means height of 0
            return 0
        else:
            return root.height  # return instance var height

    def insert(self, key):
        self.root = self._insert(self.root, key, self.count)  # Returns root of resulting tree after insertion - update it
        self.count += 1
        return self.count - 1

    def _insert(self, root: AVLTreeNode, key, val=None) -> AVLTreeNode:
        if not root:
            return AVLTreeNode(key, val, bf=0)  # If empty root this is the root of new tree
        if key < root.key:
            left_sub_root = self._insert(root.left, key, val)  # insert and update left subroot
            root.left = left_sub_root
            left_sub_root.parent = root  # assign the parent
        elif key > root.key:
            right_sub_root = self._insert(root.right, key, val)  # insert and update right subroot
            root.right = right_sub_root
            right_sub_root.parent = root
        else:
            return root  # no duplicate keys allowed; no insertion, return current root as is

        # finally, update heights and bf's of current root after insertion completed (postorder processing)
        root.height = max(self._get_height(root.left), self._get_height(root.right)) + 1
        root.bf = self._get_height(root.left) - self._get_height(root.right)
        return self.rebalance(root)  # RE-BALANCE CURRENT ROOT (if required)

    def rebalance(self, root: AVLTreeNode) -> AVLTreeNode:
        if root.bf == 2:
            if root.left.bf < 0:  # L-R
                root.left = self.rotate_left(root.left)
                return self.rotate_right(root)
            else:  # L-L
                return self.rotate_right(root)
        elif root.bf == -2:
            if root.right.bf > 0:  # R-L
                root.right = self.rotate_right(root.right)
                return self.rotate_left(root)
            else:  # R-R
                return self.rotate_left(root)
        else:
            return root  # no need to rebalance

    def rotate_right(self, root: AVLTreeNode) -> AVLTreeNode:
        pivot = root.left  # set up pointers
        tmp = pivot.right
        # 1st Move: reassign pivot's right child to root and update parent pointers
        pivot.right = root
        pivot.parent = root.parent  # pivot's parent now root's parent
        root.parent = pivot  # root's parent now pivot
        # 2nd Move: use saved right child of pivot and assign it to root's left and update its parent
        root.left = tmp
        if tmp:  # tmp can be null
            tmp.parent = root

        # update pivot's parent (manually check which one matches the root that was passed)
        if pivot.parent:
            if pivot.parent.left == root:  # if the parent's left subtree is the one to be updated
                pivot.parent.left = pivot  # assign the pivot as the new child
            else:
                pivot.parent.right = pivot  # vice-versa for right child

        # update heights and bf's using tracked heights
        root.height = max(self._get_height(root.left), self._get_height(root.right)) + 1
        root.bf = self._get_height(root.left) - self._get_height(root.right)
        pivot.height = max(self._get_height(pivot.left), self._get_height(pivot.right)) + 1
        pivot.bf = self._get_height(pivot.left) - self._get_height(pivot.right)
        return pivot  # return root of new tree

    def rotate_left(self, root: AVLTreeNode) -> AVLTreeNode:
        pivot = root.right  # set up pointers
        tmp = pivot.left
        # 1st Move: reassign pivot's right child to root and update parent pointers
        pivot.left = root
        pivot.parent = root.parent  # pivot's parent now root's parent
        root.parent = pivot  # root's parent now pivot
        # 2nd Move: use saved right child of pivot and assign it to root's left and update its parent
        root.right = tmp
        if tmp:  # tmp can be null
            tmp.parent = root

        # update pivot's parent (manually check which one matches the root that was passed)
        if pivot.parent:
            if pivot.parent.left == root:  # if the parent's left subtree is the one to be updated
                pivot.parent.left = pivot  # assign the pivot as the new child
            else:
                pivot.parent.right = pivot  # vice-versa for right child

        # update heights and bf's using tracked heights
        root.height = max(self._get_height(root.left), self._get_height(root.right)) + 1
        root.bf = self._get_height(root.left) - self._get_height(root.right)
        pivot.height = max(self._get_height(pivot.left), self._get_height(pivot.right)) + 1
        pivot.bf = self._get_height(pivot.left) - self._get_height(pivot.right)
        return pivot  # return root of new tree

    def search(self, k):
        (_upper, _lower) = (None, None)
        cur_node = self.root
        while cur_node is not None:
            if k < cur_node.key:
                _upper = cur_node
                cur_node = cur_node.left
            elif k > cur_node.key:
                _lower = cur_node
                cur_node = cur_node.right
            else:
                _upper = cur_node
                _lower = cur_node
                cur_node = None
        return _upper, _lower

    def _print_tree(self, root: AVLTreeNode):
        # (Value, idx, size)
        global tabs
        if root is None:
            print('(None)')
        else:
            print((root.key, root.val, root.height))
            tabs += 1
            print('--' * tabs, end='')
            self._print_tree(root.left)
            print('--' * tabs, end='')
            self._print_tree(root.right)
            tabs -= 1

    def print_tree(self):
        global tabs
        tabs = 0
        self._print_tree(self.root)


if __name__ == '__main__':
    seq = [1, 5, 7, -2, 10, 30]

    t = AVLTree()

    for elm in seq:
        t.insert(elm)
        t.print_tree()
        print('======')

    _upper, _lower = t.search(9)
    print(_upper, _lower)

    _upper, _lower = t.search(7)
    print(_upper, _lower)

    _upper, _lower = t.search(-5)
    print(_upper, _lower)


