import numpy as np

class Mutation:
    def __init__(self, chrom=None, allele=None, homolog_id=None, driver=True, gen=None):
        self.chrom = chrom
        self.allele = allele
        self.homolog_id = homolog_id
        self.driver = driver
        self.gen = gen
    
    def __str__(self):
        ret_str = ''
        for attr, value in vars(self).items():
            ret_str += f'{attr}: {value}, '
        return ret_str[:-2]

class CNA(Mutation):
    def __init__(self, category=None, spawn_id=None, arm=None, start=None, end=None, copies=None, gene_ids=None, **kwargs):
        self.category = category
        self.spawn_id = spawn_id
        self.arm = arm
        self.start = start
        self.gene_ids = gene_ids
        self.end = end
        self.copies = copies
        super().__init__(**kwargs)

class SNV(Mutation):
    '''
    bp: 0, 1, or 2, indicating the index of the 3 alternate base pairs ([A,C,G,T]/ref)
    '''
    def __init__(self, index=None, region=None, position=None, bp=None, **kwargs):
        self.index = index
        self.region = region
        self.position = position
        self.bp=bp
        super().__init__(**kwargs)


def derive_focal_rate(total_rate, arm_rate, chrom_rate):
    return 1 - ((1 - total_rate)/((1-arm_rate)*(1-chrom_rate)))

def draw_num_CNA_events(ncells, focal_rate, arm_rate, chromosomal_rate, WGD_rate):
    '''
    Given focal event rate, arm, chromosomal, and WGD rates, determine the number of events across all n cells of a clone
    '''
    tot_mut_rate = 1 - (1 - focal_rate) * (1 - arm_rate) * (1 - chromosomal_rate) * (1 - WGD_rate)
    nmuts = np.random.binomial(ncells, tot_mut_rate)
    return nmuts

def draw_num_SNV_events(ncells, SNV_rate):
    if SNV_rate == 0:
        return 0
    nmuts = np.random.binomial(ncells, SNV_rate)
    return nmuts

def draw_passenger_counts(cell, focal_pass_rate, SNV_pass_rate):
    num_gen = cell.length
    #npass_CNA = max(np.random.binomial(num_gen, focal_pass_rate), 1)
    npass_CNA = np.random.binomial(num_gen, focal_pass_rate)
    npass_SNV = np.random.binomial(num_gen, SNV_pass_rate)
    return npass_CNA, npass_SNV

def distribute_events(ncells, nevents_CN, nevents_SNV):
    '''
    Given n events, distribute them randomly over ncells and return an array of cell event counts, filtered if the count is > 0
    '''
    CN_counts = np.bincount(np.random.choice(ncells, nevents_CN), minlength=ncells)
    SNV_counts = np.bincount(np.random.choice(ncells, nevents_SNV), minlength=ncells)
    counts = np.stack((CN_counts, SNV_counts), axis=1)
    counts = counts[np.sum(counts, axis=1) > 0]
    return counts

def select_CNA_events(clone, nevents, focal_rate, arm_rate, chromosomal_rate, WGD_rate, focal_gain_rate, chrom_dup_rate, length_mean, mag_mean, driver=True):
    '''
    Given there are nevents, choose which events based on the relative probabiltiies
    '''
    if nevents > 0:
        rates = np.array([focal_rate, arm_rate, chromosomal_rate, WGD_rate])
        n_focal, n_arm, n_chromosomal, n_WGD = np.random.multinomial(nevents, rates / sum(rates))
        if n_WGD:
            clone.WGD()
        for i in range(n_chromosomal):
            gen_chromosomal_event(clone, chrom_dup_rate)
        for i in range(n_arm):
            gen_arm_event(clone, chrom_dup_rate)
        for i in range(n_focal):
            gen_focal_event(clone, focal_gain_rate, length_mean, mag_mean, driver=driver)

def gen_focal_event(clone, focal_gain_rate, length_mean, mag_mean, driver=True):
    '''
    Generates a focal event by drawing a chromosome and length, selecting an appropriate start location depending on driver/passenger, then selecting amp vs deletion
    '''
    poss_starts = []
    length = max(round(np.random.exponential(length_mean)), 1)

    while length > 0:
        exclude = []
        while len(poss_starts) == 0:
            chromosome = clone.genome.get_random_chrom(exclude=exclude)
            if chromosome is None: # No available chromosome
                break
            fitted_length = min(length, len(chromosome))
            if driver:
                poss_starts = clone.get_driver_start_regions(chromosome, fitted_length)
            else:
                poss_starts = clone.get_passenger_start_regions(chromosome, fitted_length)
            exclude.append(chromosome)
        if len(poss_starts) > 0:
            break
        length -= 1
        
    if len(poss_starts) == 0:
        return

    start = np.random.choice(poss_starts)
    end = start + fitted_length
    if np.random.binomial(1, focal_gain_rate):
        mag = np.random.geometric(1 / mag_mean)
        clone.focal_amplification(chromosome, start, end, mag)
    else:
        clone.focal_deletion(chromosome, start, end)

def gen_arm_event(clone, chrom_dup_rate):
    x = clone.genome.get_random_intact_arm()
    if x is not None:
        chromosome, arm = x
        if np.random.binomial(1, chrom_dup_rate):
            clone.arm_duplication(chromosome, arm)
        else:
            clone.arm_deletion(chromosome, arm)

def gen_chromosomal_event(clone, chrom_dup_rate):
    chromosome = clone.genome.get_random_intact_chrom()
    if chromosome is not None:
        if np.random.binomial(1, chrom_dup_rate):
            clone.chromosomal_duplication(chromosome)
        else:
            clone.chromosomal_deletion(chromosome)

def select_SNV_events(clone, nevents, driver=True):
    """
    Determines a set of regions to assign independent SNVs.
    """
    exclude = []
    for _ in range(nevents):
        # Start by selecting chromosome that contains appropriate regions
        poss_starts = []
        while len(poss_starts) == 0:
            chromosome = clone.genome.get_random_chrom(exclude=exclude)
            if chromosome is None: # No available chromosome
                break
            if driver:
                poss_starts = clone.get_driver_start_regions(chromosome, 1)
            else:
                poss_starts = clone.get_passenger_start_regions(chromosome, 1)
            if len(poss_starts) == 0:
                exclude.append(chromosome)

        if chromosome is None or len(poss_starts) == 0:
            return
        
        # Now choose the region and generate the event
        start = np.random.choice(poss_starts)
        gen_SNV_event(clone, chromosome, start)

def gen_SNV_event(clone, chromosome, start):
    pos = np.random.randint(clone.library.region_len)
    bp = np.random.randint(3)
    clone.add_SNV(chromosome, start, pos, bp)