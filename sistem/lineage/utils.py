import string
from itertools import product

from sistem.lineage.cell import Cell

def get_ancestors(current):
    """
    Get all clones ancestral to the current clones.

    Args:
        current (defaultdict): The current/observed clones at the time of sampling as keys with the number of sampled cells as values.

    Returns:
        list: All ancestral clones.
    """
    ancestors = set()
    for clone in current:
        par = clone.parent
        while par != None:
            ancestors.add(par)
            par = par.parent
    ancestors = list(ancestors)
    return ancestors

# Assigns names to each clone in the tree with site-specific letter code and numerical clone id
def name_clone_tree(tree, nsites):
    """
    Adds names to each clone in the tree with site-specific letter code and numerical clone id.

    Args:
        tree (Tree): The clone tree.
        nsites (int): The number of anatomical sites.

    """
    chars = list(string.ascii_uppercase)
    r = 2
    while len(chars) < nsites:
        chars = [''.join(x) for x in product(chars, repeat=r)]
        r += 1
    site_counts = {i: 0 for i in range(nsites)}
    for node in tree.iter_descendants():
        if 'diploid' not in node.name:
            node.name = chars[node.site] + '-' + str(site_counts[node.site])
            site_counts[node.site] += 1

# Assigns names to new nodes in cell tree based on those already present in clone tree
def name_cell_tree(tree):
    name_dir = {}
    for node in tree.iter_descendants():
        if node.name is None:
            par_name = node.parent.name
            if '_' in par_name:
                key = par_name.split('_')[0]
            else:
                key = par_name
            if key not in name_dir:
                name_dir[key] = 1
            node.name = key + '_' + str(name_dir[key])
            name_dir[key] += 1

def add_dummy_diploids(tree, term_gen, num_diploid=1):
    n = Cell(name='diploid', birth_gen=term_gen-1, library=tree.root.library)
    tree.root.add_child(n)
    n.inherit()
    if num_diploid > 1:
        for i in range(num_diploid):
            new_n = Cell(name=f'diploid_{i+1}', birth_gen=term_gen, library=tree.root.library)
            n.add_child(new_n)
            new_n.inherit()
