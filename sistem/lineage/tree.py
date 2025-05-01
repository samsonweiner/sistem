from collections import deque
import os
from typing import Optional

class Node:
    def __init__(
        self, 
        name=None, 
        parent=None, 
        edge_length=0, 
        **kwargs
    ):
        self.name = name
        self.parent = parent
        self.length = edge_length
        self.children = []
        self.info = {}
        #super().__init__(**kwargs)
    
    def __str__(self):
        if self.name:
            return str(self.name)
        else:
            return str(id(self))

    def __len__(self):
        return len(self.iter_leaves())
    
    #format codes --> 0: leaf names only, 1: leaf names + lengths only, 2: leaf and internal names, 3: leaf and internal names + lengths
    def _write_newick(self, terminate=True, format=3):
        if self.is_leaf():
            if format == 0 or format == 2:
                return str(self)
            elif format == 1 or format == 3:
                return str(self) + ':' + str(self.length)
        else:
            newick_str = '('
            for child in self.children:
                newick_str += child._write_newick(terminate=False, format=format) + ','
            newick_str = newick_str[:-1]
            newick_str += ')'
            if format == 2 or format == 3 or terminate:
                newick_str += str(self)
            if (format == 1 or format == 3) and not terminate :
                newick_str += ':' + str(self.length)
        if terminate:
            return newick_str + ';'
        else:
            return newick_str
    
    def is_leaf(self):
        return len(self.children) == 0
    
    def is_root(self):
        return self.parent is None

    def get_root(self):
        root = self
        while root.parent is not None:
            root = root.parent
        return root
    
    def add_child(self, node):
        self.children.append(node)
        node.parent = self
    
    def detach(self):
        if not self.is_root():
            self.parent.children.remove(self)
            self.parent = None

    #level order traversal, i.e. breadth first search from the root. Does include the root itself.
    def iter_descendants(self, include_root=True):
        nodes = deque()
        if include_root:
            nodes.append(self)
        else:
            nodes.extend(self.children)

        while nodes:
            node = nodes.popleft()
            nodes.extend(node.children)
            yield node
    
    def iter_preorder(self):
        visit_queue = deque()
        visit_queue.append(self)

        while visit_queue:
            node = visit_queue.pop()
            visit_queue.extend(node.children)
            yield node
    
    def iter_postorder(self):
        visit_queue = deque()
        return_queue = deque()
        visit_queue.append(self)

        while visit_queue:
            node = visit_queue.pop()
            return_queue.append(node)
            if not node.is_leaf():
                visit_queue.extend(node.children)
        
        while return_queue:
            node = return_queue.pop()
            yield node
    
    def iter_leaves(self):
        return [n for n in self.iter_postorder() if n.is_leaf()]

    def get_height(self):
        if self.is_leaf():
            return 0
        else:
            return max([1 + child.get_height() for child in self.children])

    def get_total_branchlen(self):
        return sum([node.length for node in self.iter_descendants()])
    
    def is_ancestor(self, other):
        for node in self.iter_descendants():
            if node == other:
                return True
        return False
    
class Tree:
    def __init__(
        self, 
        root: Optional[Node] = None, 
        newick: Optional[str] = None
    ):
        self.root = root
        if newick:
            if os.path.exists(newick):
                with open(newick) as f:
                    newick_str = f.readline().strip()
                    self.root = str_to_newick(newick_str)
            else:
                self.root = str_to_newick(newick)
    
    
    def __repr__(self):
        return self.root._write_newick(format=3)

    def __len__(self):
        return len([leaf for leaf in self.iter_leaves()])

    def __iter__(self):
        return self.iter_leaves()

    def print_newick(self, format=3):
        newick_str = self.root._write_newick(format=format)
        print(newick_str)
    
    def write_newick(self, format=3):
        newick_str = self.root._write_newick(format=format)
        return newick_str

    #post order traversal, i.e. left,right,middle.
    def iter_postorder(self):
        return self.root.iter_postorder()
    
    #level order traversal, i.e. breadth first search from the root. Does include the root itself.
    def iter_descendants(self, include_root=True):
        return self.root.iter_descendants(include_root=include_root)
    
    def iter_preorder(self):
        return self.root.iter_preorder()
    
    def traverse(self):
        return self.iter_descendants()
    
    def iter_leaves(self):
        for node in self.iter_descendants():
            if node.is_leaf():
                yield node
    
    def iter_ancestors(self):
        for node in self.iter_descendants():
            if not node.is_leaf():
                yield node
    
    def iter_path_to_node(self, node):
        path = []
        cur_node = node
        while not cur_node.is_root():
            path.append(cur_node)
            cur_node = cur_node.parent
        path.append(self.root)
        path.reverse()
        return path

    def has_leaf_names(self):
        for leaf in self.iter_leaves():
            if leaf.name is None:
                return False
        return True

    def set_leaf_names(self):
        count = 1
        for leaf in self.iter_leaves():
            if leaf.name is None:
                leaf.name = 'cell' + str(count)
                count += 1
    
    def set_node_names(self):
        count = 0
        for n in self.iter_postorder():
            if n.name == None:
                n.name = f'A{count}'
                count += 1

    def get_tree_height(self):
        return self.root.get_height()

    def get_total_branchlen(self):
        return self.root.get_total_branchlen()

    def find(self, node_name):
        for node in self.iter_descendants():
            if node.name == node_name:
                return node
        return None

    def get_mutations(self):
        mutations = []
        for node in self.iter_descendants():
            mutations += node.events
        return mutations

    def save(self, file_path, format=3):
        newick_str = self.root._write_newick(format=format)
        f = open(file_path, 'w+')
        f.write(newick_str)
        f.close()

    def copy_structure(self):
        ts = self._write_newick(format=3)
        new_t = Tree(newick=ts)
        return new_t
    

    # Function which resolves all unifurcating nodes in a tree by merging with their children. Updates branch lengths.
    def resolve_unifurcations(self):
        nodes = [self.root]
        while nodes:
            node = nodes.pop(0)
            children = node.children
            if len(children) == 1:
                child = children[0]
                child.detach()
                node.children = child.children
                for gc in node.children:
                    gc.parent = node

                child.children = []
                node.length += child.length
                nodes.append(node)

            else:
                for c in node.children:
                    nodes.append(c)
    
    def set_branchlengths(self):
        """
        Sets the branch lengths according to difference in birth gens.
        """
        for node in self.iter_postorder():
            if not node.is_root():
                cur_gen = node.birth_gen
                par_gen = node.parent.birth_gen
                node.length = cur_gen - par_gen

def str_to_newick(newick_str):
    split_str = newick_str[:-1].split(',')

    cur_node = None
    nodes = []

    for chunk in split_str:
        while chunk[0] == '(':
            new_node = Node()
            if cur_node:
                cur_node.add_child(new_node)
            cur_node = new_node
            chunk = chunk[1:]
        rest = chunk.split(')')
        if ':' in rest[0]:
            idx = rest[0].index(':')
            try:
                edge_len = float(rest[0][idx+1:])
            except:
                edge_len = 0
            child_node = Node(name=rest[0][:idx], edge_length=edge_len)
            cur_node.add_child(child_node)
        else:
            child_node = Node(name=rest[0])
            cur_node.add_child(child_node)
        if len(rest) > 1:
            for part in rest[1:]:
                if ':' in part:
                    idx = part.index(':')
                    cur_node.name = part[:idx]
                    try:
                        edge_len=float(part[idx+1:])
                    except:
                        edge_len=0
                    cur_node.length = edge_len
                else:
                    cur_node.name = part
                    cur_node.length = 0
                if not cur_node.is_root():
                    cur_node = cur_node.parent
    return cur_node


