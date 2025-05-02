import copy

from sistem.selection import Attributes
from sistem.selection.arm_library import BaseArmLibrary
from sistem.genome.utils import ConvertToSNVGenome
from sistem.lineage.mutation import CNA, SNV
import sistem.lineage.mutation as mut
from sistem.utilities.utilities import combine_dicts, get_reg_id

def evolve_tree_with_passengers(tree, igenome, focal_pass_rate, SNV_pass_rate, focal_gain_rate, length_mean, mag_mean, no_leaves=False):
    """
    Master function for re-evolving a tree with passenger mutations. Tree must be resolved with appropriate branch lengths.

    Args:
        tree (Tree): A tree object made of Cell objects.
        igenome (Genome): Initial starting genome
        focal_gain_rate (float): The probability of a passenger CNA being an amplification (deletion rate is 1-focal_gain_rate).
        length_mean (int): The mean length of passenger CNAs in terms of number of regions.
        mag_mean (float): The mean number of additional copies gained during amplifications.
        no_leaves (Bool): If toggled, does not draw passenger mutations for leaf nodes
    """
    if tree.root.library.is_driver_SNV_model:
        tree.root.genome = copy.deepcopy(igenome)
    else:
        altgenome = ConvertToSNVGenome(igenome)
        tree.root.genome = copy.deepcopy(altgenome)

    reevolve_tree(tree)

    if not tree.root.library.is_driver_region_model and isinstance(tree.root.library, BaseArmLibrary):
        get_blacklisted_regions(tree.root)

    for cell in tree.iter_descendants(include_root=False):
        if no_leaves and cell.is_leaf():
            cell.inherit()
        else:
            inherit_with_passengers(cell)
            npass_CNA, npass_SNV = mut.draw_passenger_counts(cell, focal_pass_rate, SNV_pass_rate)
            for _ in range(npass_CNA):
                mut.gen_focal_event(cell, focal_gain_rate, length_mean, mag_mean, driver=False)
            if npass_SNV > 0:
                mut.select_SNV_events(cell, npass_SNV, driver=False)

def reevolve_tree(tree):
    """
    Re-applies driver mutations across tree, recording mutated subsequence of focal events along the way

    Args:
        tree (Tree): A tree object made of Cell objects.
    """
    for node in tree.iter_descendants():
        if not node.is_root():
            node.inherit()
        for event in node.events:
            apply_event(node, event)

def apply_event(cell, event):
    if isinstance(event, CNA):
        if event.category == 'WGD':
            cell.WGD(record_event=False)
        else:
            cur_chromosome = cell.genome.find(event.chrom, event.homolog_id)
            if event.category == 'focal_amplification':
                event.gene_ids = cur_chromosome.seq[event.start:event.end]
                cell.focal_amplification(cur_chromosome, event.start, event.end, event.copies, record_event=False)
            elif event.category == 'focal_deletion':
                event.gene_ids = cur_chromosome.seq[event.start:event.end]
                cell.focal_deletion(cur_chromosome, event.start, event.end, record_event=False)
            elif event.category == 'arm_duplication':
                cell.arm_duplication(cur_chromosome, event.arm, record_event=False)
            elif event.category == 'arm_deletion':
                cell.arm_deletion(cur_chromosome, event.arm, record_event=False)
            elif event.category == 'chromosomal_duplication':
                cell.chromosomal_duplication(cur_chromosome, record_event=False)
            elif event.category == 'chromosomal_deletion':
                cell.chromosomal_deletion(cur_chromosome, record_event=False)
    elif isinstance(event, SNV):
        cur_chromosome = cell.genome.find(event.chrom, event.homolog_id)
        event.region = cur_chromosome.seq[event.index]
        cell.add_SNV(cur_chromosome, event.index, event.position, event.bp, record_event=False)

def get_blacklisted_regions(node):
    """
    Function for arm library model. Finds regions that are mutated with focal driver events in the agent-based model to ensure that passenger mutations do not interfere with these events. Sets as a 'blacklisted_regions' attribute.

    Args:
        node (Node): Clone/Cell/Node object
    
    Returns:
        blacklisted_regions (dict): Dictionary of blacklisted regions.
    """
    child_regions = []
    for child in node.children:
        child_regions.append(get_blacklisted_regions(child))
    blacklisted_regions = combine_dicts(*child_regions)

    for event in node.events:
        if isinstance(event, CNA):
            if event.category == 'focal_amplification' or event.category == 'focal_deletion':
                if event.chrom not in blacklisted_regions:
                    blacklisted_regions[event.chrom] = {}
                if event.homolog_id not in blacklisted_regions[event.chrom]:
                    blacklisted_regions[event.chrom][event.homolog_id] = set()
                blacklisted_regions[event.chrom][event.homolog_id].update([get_reg_id(r) for r in event.gene_ids])
            elif event.category == 'chromosomal_duplication' or event.category == 'arm_duplication':
                if event.chrom in blacklisted_regions:
                    if event.spawn_id in blacklisted_regions[event.chrom]:
                        if event.homolog_id in blacklisted_regions[event.chrom]:
                            blacklisted_regions[event.chrom][event.homolog_id].update(blacklisted_regions[event.chrom][event.spawn_id])
                        else:
                            blacklisted_regions[event.chrom][event.homolog_id] = blacklisted_regions[event.chrom][event.spawn_id].copy()
            elif event.category == 'WGD':
                for chrname in node.genome.chrom_names:
                    if chrname in blacklisted_regions:
                        n = int(len(node.genome[chrname])/2)
                        for spawn_id in range(n, len(node.genome[chrname])):
                            if spawn_id in blacklisted_regions[chrname]:
                                homolog_id = spawn_id - n
                                if homolog_id in blacklisted_regions[chrname]:
                                    blacklisted_regions[chrname][homolog_id].update(blacklisted_regions[chrname][spawn_id])
                                else:
                                    blacklisted_regions[chrname][homolog_id] = blacklisted_regions[chrname][spawn_id].copy()

    node.attributes = Attributes(blacklisted_regions=blacklisted_regions)
    return blacklisted_regions


def inherit_with_passengers(cell):
    """
    Inherits genome from parent while correcting driver mutations indices based on new passenger mutations.

    Args:
        cell (Cell): The cell to inherit.
    """
    assert cell.parent is not None and cell.parent.genome is not None, "Cannot inherit if parent and/or parent genome non-existent"
    cell.inherit()
    # If there are driver events, there is always only one
    for event in cell.events:
        if isinstance(event, CNA):
            if event.category == 'WGD':
                cell.WGD(record_event=False)
            else:
                cur_chromosome = cell.genome.find(event.chrom, event.homolog_id)
                if event.category == 'chromosomal_duplication':
                    cell.chromosomal_duplication(cur_chromosome, record_event=False)
                elif event.category == 'arm_duplication':
                    cell.arm_duplication(cur_chromosome, event.arm, record_event=False)
                elif event.category == 'chromosomal_deletion':
                    cell.chromosomal_deletion(cur_chromosome, record_event=False)
                elif event.category == 'arm_deletion':
                    cell.arm_deletion(cur_chromosome, event.arm, record_event=False)
                elif event.category == 'focal_amplification' or event.category == 'focal_deletion':
                    inherit_with_focal(event, cur_chromosome, cell.library)
        elif isinstance(event, SNV):
            cur_chromosome = cell.genome.find(event.chrom, event.homolog_id)
            inherit_with_SNV(event, cur_chromosome, cell.library)

def inherit_with_focal(event, chromosome, library):
    """
    Inherits chromosome from parent and applies driver event with updated region indices. 
    """

    # Logic for if the library uses driver regions
    if library.is_driver_region_model:
        # get number of regions that are deleted before and after the first/last driver region in mutated segment
        s1, r1, s2, r2 = get_first_last_driver(event.gene_ids, library, event.chrom)
        
        # find regions in parent chromosome, make sure regions on either side are still passengers
        #nchromosome = source_chromosome.copy()
        i1 = chromosome.seq.index(r1)
        i2 = chromosome.seq.index(r2) + 1
        for _ in range(s1):
            if i1-1 == 0 or chromosome.seq[i1 - 1] in library.drivers[event.chrom]:
                break
            i1 -= 1
        for _ in range(s2):
            if i2 >= len(chromosome.seq) or chromosome.seq[i2 + 1] in library.drivers[event.chrom]:
                break
            i2 += 1
    
    #if the library does not use driver regions (e.g. arm model), then we preserve the entire mutated subsequence
    else:
        r1,r2 = event.gene_ids[0], event.gene_ids[-1]
        i1 = chromosome.seq.index(r1)
        i2 = chromosome.seq.index(r2) + 1

    # update event with new indices
    event.start = i1
    event.end = i2
    event.gene_ids = chromosome.seq[event.start:event.end]
    if event.category == 'focal_amplification':
        chromosome.amplify(i1, i2, event.copies)
    elif event.category == 'focal_deletion':
        chromosome.delete(i1, i2)

def get_first_last_driver(segment, library, chrname):
    """
    Get the number of continous passenger regions at the beginning and end of a segment. Only used in the driver gene model.

    Args:
        segment (list): A subset of a chromosome.seq corresponding to the mutated segment.
        library (Library): Cell Library.
        chrname (str): Name of the mutated chromosome.

    Returns:
        s1 (int): Number of regions before the first driver region in the mutated sequence.
        r1 (int): Region id at s1.
        s2 (int): Number of regions before the last driver region in the mutated sequence.
        r2 (int): Region id at s2.
    """
    for i,r in enumerate(segment):
        if get_reg_id(r) in library.drivers[chrname] + library.essential[chrname]:
            s1 = i
            r1 = r
            break
    #for i,r in reversed(zip(enumerate(segment))):
    for i,r in enumerate(reversed(segment)):
        if get_reg_id(r) in library.drivers[chrname] + library.essential[chrname]:
            s2 = i
            r2 = r
            break
    return s1, r1, s2, r2


def inherit_with_SNV(event, chromosome, library):
    if library.is_driver_region_model and library.is_driver_SNV_model:
        index = chromosome.seq.index(event.region)
        event.index = index
        chromosome.add_SNV(index, event.position, event.bp)