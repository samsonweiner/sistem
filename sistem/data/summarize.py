import numpy as np

def extract_largest_clones(gs):
    '''
    Finds the largest clone in each anatomical site based on popsize and returns it
    '''
    largest = []
    for s,clones in gs.clones.items():
        if len(clones):
            c = max(clones, key = lambda x: x.popsize)
            largest.append(c)
        else:
            largest.append(None)
    return largest

def arm_CNratios_single(cell):
    CNratios = {}
    for chrname in cell.genome.chrom_names:
        CNratios[chrname] = []
        for a in [0, 1]:
            CNrat = cell.attributes.arm_counts[chrname][a]/(cell.library.arm_sizes[chrname][a])
            CNratios[chrname].append(CNrat)
    return CNratios

def arm_CNratios_bulk(clones):
    popsizes = [clone.popsize for clone in clones]
    library = clones[0].library
    CNratios = {}
    for chrname in clones[0].genome.chrom_names:
        CNratios[chrname] = []
        for a in [0, 1]:
            arm_counts = [clone.attributes.arm_counts[chrname][a] for clone in clones]
            avg_arm_count = np.average(arm_counts, weights=popsizes)
            CNrat = avg_arm_count / library.arm_sizes[chrname][a]
            CNratios[chrname].append(CNrat)
    return CNratios