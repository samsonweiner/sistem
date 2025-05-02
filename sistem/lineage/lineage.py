import numpy as np

from sistem.lineage.cell import Cell
from sistem.lineage.tree import Tree
from sistem.lineage.utils import name_cell_tree

def create_clone_tree(observed):
    """
    Creates a tree from observed clones tracing backwards in time. Root will always be the starting clone.

    Args:
        observed (defaultdict): The observed clones at the time of sampling as keys with the number of sampled cells as values.
        nsites (int): The number of anatomical sites.

    Returns:
        Tree: An unresolved clone tree.
    """
    queue = list(observed.keys())
    nodes = [clone for clone in queue]
    serviced = set()
    while queue:
        cur = queue.pop(0)
        par = cur.parent
        serviced.add(cur)
        if par != None:
            if par not in nodes:
                nodes.append(par)
            par.children.append(cur)
            if par not in serviced and par not in queue:
                queue.append(par)

    for clone in nodes:
        if clone.parent == None:
            root = clone
    t = Tree(root=root)
    #name_clone_tree(t, nsites)
    t.set_branchlengths()
    return t

def create_converted_tree(tree, observed):
    clones = [clone for clone in tree.iter_postorder()]
    converter = {}
    for clone in clones:
        n = Cell(library=clone.library, site=clone.site, birth_gen=clone.birth_gen, name=clone.name)
        n.birth_counts = clone.birth_counts
        n.events = clone.events
        converter[clone] = n
    for clone in clones:
        cell = converter[clone]
        for child in clone.children:
            cell.children.append(converter[child])
        if clone.parent is not None:
            cell.parent = converter[clone.parent]

    new_tree = Tree(root=converter[tree.root])
    new_tree.set_branchlengths()
    observed_converted = {converter[k]: v for k,v in observed.items()}
    return new_tree, observed_converted

# Creates an individual Cell object for each observed cell of a clone and adds them as children
def create_singlecell_tree(tree, observed, term_gen):
    clones = [clone for clone in tree.iter_postorder()]
    converter = {}
    for clone in clones:
        n = Cell(library=clone.library, site=clone.site, birth_gen=clone.birth_gen, name=clone.name)
        n.birth_counts = clone.birth_counts
        n.events = clone.events
        converter[clone] = n
    for clone in clones:
        cell = converter[clone]
        for child in clone.children:
            cell.children.append(converter[child])
        if clone.parent is not None:
            cell.parent = converter[clone.parent]

    new_tree = Tree(root=converter[tree.root])
    new_tree.set_branchlengths()
    observed_converted = {converter[k]: v for k,v in observed.items()}

    for cell in new_tree.iter_descendants():
        if cell in observed_converted:
            nsampled = observed_converted[cell]
            for _ in range(nsampled):
                n = Cell(library=cell.library, site=cell.site, birth_gen=term_gen)
                cell.children.append(n)
                n.parent = cell
                
    return new_tree

# Takes as input a clone tree, and creates node objects for each individual sampled cell and resolves the topology
def resolve_multifurcation(tree, name_tree=True):
    nodes = [node for node in tree.iter_descendants() if not node.is_leaf()]
    while nodes:
        clone = nodes.pop(0)
        if len(clone.children) > 1:
            resolve_multifurcation_node(clone)
    if name_tree:
        name_cell_tree(tree)
    tree.set_branchlengths()

# Function which resolves multifurcating nodes by randomly selecting replication groups from the larger population. 
def resolve_multifurcation_node(clone, verbose=False):
    children = [c for c in clone.children]
    assert len(children) > 1, "Cannot resolve multifurcation with fewer than 2 children"

    groups = []
    t0 = clone.birth_gen
    t = max([clone.birth_gen for clone in children])

    if verbose:
        print(f'Resolving {clone.name}, {clone.birth_gen} birth_gen, {len(clone.children)} children')
    #diverge_count = {tt: 0 for tt in clone.birth_counts.keys()}
    #for child in children:
    #    diverge_count[child.birth_gen] += 1

    while t > t0:
        # Draw cell from greater pop for existing replication groups and add nodes if coalescence occurs
        # Add any children born in gen t to the set of replication groups
        for i in list(range(len(children)))[::-1]:
            if children[i].birth_gen == t:
                groups.append(children[i])
                children.pop(i)

        if verbose:
            print(f'{t}: {clone.birth_counts[t]}, {len(groups)}')
                
        assert len(groups) <= 2*clone.birth_counts[t], 'Unable to coalesce: not enough living cells from previous generation'
        pool = np.repeat(np.arange(clone.birth_counts[t]), 2)
        np.random.shuffle(pool)
        #pool = pool[diverge_count[t]:]
        assignments = pool[:len(groups)]
        draws = {}
        new_groups, remove_idxs = [], set()
        for i,b in enumerate(assignments):
            if b in draws:
                j = draws[b]
                n = Cell(library=clone.library, site=clone.site, birth_gen=t)
                n.children = [groups[i], groups[j]]
                groups[i].parent = groups[j].parent = n
                new_groups.append(n)
                remove_idxs.update((i,j))
            else:
                draws[b] = i
        
        groups = [x for i,x in enumerate(groups) if i not in remove_idxs] + new_groups

        t -= 1
    
    # Connect parent to new children in groups. If only one remaining group (creating unifurcation), leave as is for now.
    while len(groups) > 2:
        np.random.shuffle(groups)
        g1 = groups.pop(0)
        g2 = groups.pop(0)
        n = Cell(library=clone.library, site=clone.site, birth_gen=t)
        n.children = [g1, g2]
        g1.parent = g2.parent = n
        groups.append(n)

    clone.children = []
    for g in groups:
        clone.children.append(g)
        g.parent = clone
    
    if verbose:
        print('')

def merge_and_resolve_clone_tree(clone_tree, observed, nsites, term_gen, min_mut_fraction):
    """
    
    """

    # First add one child cell to each observed clone, then resolve multifurcations
    sample_counts = {}
    for clone, count in observed.items():
        n = Cell(library=clone.library, site=clone.site, birth_gen=term_gen)
        clone.add_child(n)
        sample_counts[n] = count
    resolve_multifurcation(clone_tree)

    # For each node, count the total number of samples cells in each site which are descendant from that node
    subtree_counts = {node: np.zeros(nsites) for node in clone_tree.traverse()}
    for node in clone_tree.iter_postorder():
        if node.is_leaf():
            subtree_counts[node][node.site] += sample_counts[node]
        else:
            for child in node.children:
                subtree_counts[node] += subtree_counts[child]

    # Helper function which, given a node and a site, removes all descendant clones (possibly itself) which belong to the site
    def remove_subclones_in_site(node, site, exclude=[]):
        if node.is_leaf():
            if node.site == site:
                if node not in exclude:
                    del sample_counts[node]
                    node.detach()
        else:
            children = node.children[:]
            for child in children:
                remove_subclones_in_site(child, site, exclude=exclude)
            if len(node.children) == 0:
                node.detach()
    
    # Helper function which, for each site, looks for subtrees in which the total number of observed cells from that site is less than some minimum threshold (min_cell_per_site). If so, will either prune the subtree or merge it with other clones.
    # min_cell_per_site: 
    def merge_clones(node, site, min_cell_per_site):
        if node.is_leaf():
            # If at this point, the leaf should already exceed the threshold, otherwise it would've been removed.
            return
        if subtree_counts[node][site] >= min_cell_per_site:
            if len(node.children) == 1:
                merge_clones(node.children[0], site, min_cell_per_site)
            else:
                assert len(node.children) == 2, 'Somehow more than 2 children'
                left, right = node.children[0], node.children[1]
                # case 1: both subtrees exceeds the lower bound
                if subtree_counts[left][site] >= min_cell_per_site and subtree_counts[right][site] >= min_cell_per_site:
                    merge_clones(left, site, min_cell_per_site)
                    merge_clones(right, site, min_cell_per_site)
                # case 2: only the left subtree passes, prune the right
                elif subtree_counts[left][site] >= min_cell_per_site:
                    remove_subclones_in_site(right, site)
                    merge_clones(left, site, min_cell_per_site)
                # case 3: only the right subtree passes, prune the left
                elif subtree_counts[right][site] >= min_cell_per_site:
                    remove_subclones_in_site(left, site)
                    merge_clones(right, site, min_cell_per_site)
                # case 4: neither passes
                else:
                    # Need to combine in a single clone. For now, we'll just pretend all of the other clones belong to 'largest'.
                    assert subtree_counts[left][site] + subtree_counts[right][site] >= min_cell_per_site, 'This should be the only remaining case'
                    leaves = [leaf for leaf in node.iter_leaves() if leaf.site == site]
                    leaves.sort(key = lambda x: sample_counts[x], reverse=True)
                    largest = leaves[0]
                    sample_counts[largest] = sum([sample_counts[leaf] for leaf in leaves])
                    #for leaf in leaves[1:]:
                    #    del sample_counts[leaf]
                    remove_subclones_in_site(node, site, exclude=[largest])

        # No possible way to overcome threshold, so just remove all the clones in the site in the subtree
        else:
            remove_subclones_in_site(node, site)

    # Counts the total number of observed cells in each site (equal to [ncells_prim, ncells_meta, ncells_meta, ...])
    tot_observed = [0 for _ in range(nsites)]
    for leaf,count in sample_counts.items():
        tot_observed[leaf.site] += count

    # Compute min_cell_per_site as the total number of observed cells in that site * a minimum mutation fraction.
    for s in range(nsites):
        min_cell_per_site = int(tot_observed[s]*min_mut_fraction)
        merge_clones(node, s, min_cell_per_site)

    #add_dummy_diploid(clone_tree, term_gen)
    #collapse_tree(clone_tree)
    #name_clone_tree(clone_tree, nsites)
    return sample_counts

def collapse_tree(tree):
    def collapse_tree_helper(node):
        if len(node.children) == 1:
            child = node.children[0]
            child.events = node.events + child.events
            #if not child.is_leaf():
            #    child.site = node.site
            child.detach()
            if node.parent == None:
                tree.root = child
            else:
                parent = node.parent
                node.detach()
                parent.add_child(child)
            collapse_tree_helper(child)
        elif len(node.children) > 1:
            children = node.children[:]
            for child in children:
                collapse_tree_helper(child)
    
    collapse_tree_helper(tree.root)
    tree.set_branchlengths()