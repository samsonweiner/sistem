import numpy as np
from scipy.stats import expon
from itertools import chain
from collections import defaultdict
import logging
import pickle
import copy
import os

from typing import Optional, Union

from sistem.utilities.IO import setup_logger
from sistem.genome.genome import init_diploid_genome
from sistem.lineage import Cell, Clone
import sistem.lineage.mutation as mut
from sistem.lineage.lineage import create_clone_tree, create_converted_tree, create_singlecell_tree, resolve_multifurcation, merge_and_resolve_clone_tree, collapse_tree
from sistem.lineage.reevolve import evolve_tree_with_passengers
from sistem.lineage.utils import add_dummy_diploids
from sistem.anatomy.utils import get_num_migrations, draw_migrations, generate_migration_graph, save_migration_graph
from sistem.anatomy import BaseAnatomy
from sistem.data.profiles import save_mutations_from_tree, save_CNPs, save_SNVs, site_CN_averages
from sistem.parameters import Parameters, fill_params

class GrowthSimulator:
    """The main class used to run the simulator. 

    Args:
        anatomy (Anatomy): The Anatomy object which has been initialized. Required.

    Attributes:
        clones (Dict[list]): A dictionary where keys are site numbers (0-nsites-1), and values are a list of clones currently present in the site. 
        site_counts (array): Contains the total cell count of each site.
        gen (int): The current generation of the simulator.
    """
    def __init__(self, anatomy):
        """Constructor method

        """
        self.anatomy = self._check_anatomy(anatomy)

        self.clones = {i: [] for i in range(self.anatomy.nsites)}
        self.site_counts = np.zeros(self.anatomy.nsites, dtype=int)
        self.gen = 0

        self._empty_sites = set()
        self._clone_ids = [0 for _ in range(self.anatomy.nsites)]
        self._igenome = init_diploid_genome(self.anatomy.libraries[0])
        self._init_clone = Clone(name=self._assign_clone_name(0), genome=copy.deepcopy(self._igenome), library=self.anatomy.libraries[0], popsize=self.anatomy.N0, site=0, birth_gen=0)
        self.clones[0].append(self._init_clone)
        self.site_counts[0] += self._init_clone.popsize

        self._lifespan_gens = {i: defaultdict(list) for i in range(self.anatomy.nsites)}
        self._lifespan_gens[0][1].append(self._init_clone)

        # Data structures/variables used in sampling and phylogeny construction
        self._observed = defaultdict(int)
        self._clone_tree = None

    def __iter__(self):
        return chain.from_iterable(self.clones.values())

    def _check_anatomy(self, anatomy):
        if not isinstance(anatomy, BaseAnatomy):
            raise TypeError("Argument `anatomy' must be an instance of an Anatomy class.")
        if anatomy.initialized is False:
            raise ValueError("Anatomy object has not been properly initialized. Try running anatomy.initialize_distances method.")
        return anatomy
    
    def _assign_clone_name(self, site):
        name = self.anatomy.site_ids[site] + str(self._clone_ids[site])
        self._clone_ids[site] += 1
        return name
    
    def _get_mean_fitness(self):
        for clone in self:
            clone.update_fitness()
        mean_fits = []
        for s,clones in self.clones.items():
            if len(clones) == 0:
                mean_fits.append(0)
            else:
                mean_fits.append(np.average([clone.fitness for clone in clones], weights=[clone.popsize for clone in clones]))
        return mean_fits
    
    def _create_new_clones(self, parent_clone, site, distr_events, focal_rate, arm_rate, chromosomal_rate, WGD_rate, focal_gain_rate, chrom_dup_rate, length_mean, mag_mean, lifespan_mean):
        """
        Given array distr_events from draw_num_CNA_events, creates a new clone for each element with event count equal to the element
        """
        nviolations = 0
        for x,y in distr_events:
            clone = parent_clone.replicate(name=self._assign_clone_name(site))
            if x > 0:
                mut.select_CNA_events(clone, x, focal_rate, arm_rate, chromosomal_rate, WGD_rate, focal_gain_rate, chrom_dup_rate, length_mean, mag_mean)
            if y > 0:
                mut.select_SNV_events(clone, y)
            if clone.is_viable():
                self.clones[site].append(clone)
                if lifespan_mean > 1:
                    lifespan = max(1, round(expon.rvs(scale=2)))
                else:
                    lifespan = 1
                self._lifespan_gens[site][self.gen + lifespan].append(clone)
            else:
                nviolations += 1
                clone.die()
        return nviolations
    
    # PARALLELIZE?
    def _cycle_birthdeath(self, focal_driver_rate, arm_rate, chromosomal_rate, WGD_rate, focal_gain_rate, chrom_dup_rate, length_mean, mag_mean, SNV_driver_rate, lifespan_mean):
        totbirths, totdeaths = [0 for i in range(self.anatomy.nsites)], [0 for i in range(self.anatomy.nsites)]
        mean_fits = self._get_mean_fitness()
        
        # Driver mutations
        for s,m in enumerate(mean_fits):
            #clones = self.clones[s][:]
            clones = self._lifespan_gens[s][self.gen]
            if len(clones) > 1:
                print(self.gen, clones)
            if len(clones) > 0:
                Ntot, Etot = self.site_counts[s], self.anatomy.exp_pop[s][self.gen]
                #for clone in clones:
                while len(clones) > 0:
                    clone = clones.pop(0)
                    nbirths = clone.get_num_births(m, Ntot, Etot)
                    clone.birth_counts[self.gen] = nbirths
                    ndeaths = clone.popsize - nbirths
                    if nbirths > 0:
                        clone.popsize = 2*nbirths
                        clone_focal_rate = focal_driver_rate*clone.get_ploidy()/2
                        nevents_CN = mut.draw_num_CNA_events(nbirths, clone_focal_rate, arm_rate, chromosomal_rate, WGD_rate)
                        nevents_SNV = mut.draw_num_SNV_events(nbirths, SNV_driver_rate)
                        if nevents_CN > 0 or nevents_SNV > 0:
                            distr_events = mut.distribute_events(nbirths, nevents_CN, nevents_SNV)
                            nviolations = self._create_new_clones(clone, s, distr_events, clone_focal_rate, arm_rate, chromosomal_rate, WGD_rate, focal_gain_rate, chrom_dup_rate, length_mean, mag_mean, lifespan_mean)
                            clone.popsize -= len(distr_events)
                            self.site_counts[s] -= nviolations

                        if clone.popsize <= 0:
                            clone.die()
                            self.clones[s].remove(clone)
                        else:
                            if lifespan_mean > 1:
                                lifespan = max(1, round(expon.rvs(scale=2)))
                            else:
                                lifespan = 1
                            self._lifespan_gens[s][self.gen + lifespan].append(clone)
                    else:
                        clone.die()
                        self.clones[s].remove(clone)
                    self.site_counts[s] += nbirths - ndeaths
                    totbirths[s] += nbirths
                    totdeaths[s] += ndeaths
        
        # Passenger mutations
        #for clone in self:
        #    add_passengers(clone, focal_pass_rate, SNV_pass_rate)

        return mean_fits, totbirths, totdeaths

    def _cycle_migrations(self, t_max, pattern):
        # First figure out what migrations are happening, then afterwards do the shuffling
        totmigrations = [0 for i in range(self.anatomy.nsites)]
        if self.anatomy.nsites == 1:
            return totmigrations

        multi_source = True
        if 'S' in pattern:
            multi_source = False
        
        moves = defaultdict(list)
        for s,clones in self.clones.items():
            if len(clones) == 0:
                self._empty_sites.add(s)
            for clone in clones[:]:
                nmigrations = get_num_migrations(clone, self.anatomy.get_total_transition_probability(clone))
                if nmigrations > 0:
                    migrations = draw_migrations(nmigrations, self.anatomy.get_oneway_transition_probabilities(clone))
                    for a,occurs in enumerate(migrations):
                        if occurs == 1:
                            migration_size = 1
                            if multi_source:
                                migration_size = self.anatomy.N0

                            new_clone = clone.replicate(name=self._assign_clone_name(a), site=a, popsize=migration_size, library=self.anatomy.libraries[a])
                            moves[a].append(new_clone)
                            self.site_counts[a] += migration_size
                            self.site_counts[s] -= migration_size
                            clone.popsize -= migration_size

                    if clone.popsize <= 0:
                        self.clones[s].remove(clone)
                    totmigrations[s] += sum(migrations)
        
        for s,new_clones in moves.items():
            # For newly seeded sites, need to initialize growth functions
            if s in self._empty_sites:
                fits = []
                for clone in new_clones:
                    clone.update_fitness()
                    fits.append(clone.fitness)
                mean_fit = np.mean(fits)
                self.anatomy.init_pop_dynamic(s, self.gen, t_max, cur_fit=mean_fit, library=self.anatomy.libraries[s])
                self._empty_sites.remove(s)
            self.clones[s].extend(new_clones)
        
        return totmigrations

    # TODO: implement min detectable
    def _terminate(self, t_max, min_detectable):
        if sum(self.site_counts) == 0:
            logging.info('\nTERMINATE: No remaining living cells')
            return True
        elif self.gen > t_max:
            logging.info('\nTERMINATE: Max gen reached')
            return True
        elif np.all(self.site_counts >= min_detectable):
            logging.info('\nTERMINATE: Cell populations reached detectabled sizes in all sites')
            return True
        else:
            return False
    
    def _log_current_gen(self, mean_fits, totbirths, totdeaths, totmigrations):
        logging.info(f'\nGen{self.gen}')
        logging.info(f'P\t{len(self.clones[0])}\t{self.site_counts[0]}\t{mean_fits[0]}\t{totbirths[0]}\t{totdeaths[0]}\t{totmigrations[0]}')
        for s in range(1, self.anatomy.nsites):
            logging.info(f'M{s}\t{len(self.clones[s])}\t{self.site_counts[s]}\t{mean_fits[s]}\t{totbirths[s]}\t{totdeaths[s]}\t{totmigrations[s]}')

    def simulate_agents(
        self, 
        params: Optional[Parameters] = None,
        t_max: Optional[int] = None, 
        lifespan_mean: Optional[Union[int, float]] = None,
        min_detectable: Optional[int] = None, 
        focal_driver_rate: Optional[float] = None, 
        arm_rate: Optional[float] = None, 
        chromosomal_rate: Optional[float] = None, 
        WGD_rate: Optional[float] = None, 
        focal_gain_rate: Optional[float] = None, 
        chrom_dup_rate: Optional[float] = None, 
        length_mean: Optional[Union[int, float]] = None, 
        mag_mean: Optional[Union[int, float]] = None, 
        SNV_driver_rate: Optional[float] = None, 
        log_path: Optional[str] = None,
        pattern: str = ''
    ):
        """Begins the agent-based phase of the simulator. Runs until a minimum detectable number of cells is reached in each site, the max generation is reached, or no cells are left in any site. Only simulates driver mutations. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Parameters, optional):
            t_max (int, optional):
            lifespan_mean (int, float, optional):
            min_detectable (int, optional):
            focal_driver_rate (float, optional):
            arm_rate (float, optional):
            chromosomal_rate (float, optional):
            WGD_rate (float, optional):
            focal_gain_rate (float, optional):
            chrom_dup_rate (float, optional):
            length_mean (int, float, optional):
            mag_mean (int, float, optional):
            SNV_driver_rate (float, optional):
            log_path (str, optional): 
        """
        params = fill_params(params, t_max=t_max, lifespan_mean=lifespan_mean, min_detectable=min_detectable, focal_driver_rate=focal_driver_rate, arm_rate=arm_rate, chromosomal_rate=chromosomal_rate, WGD_rate=WGD_rate, focal_gain_rate=focal_gain_rate, chrom_dup_rate=chrom_dup_rate, length_mean=length_mean, mag_mean=mag_mean, SNV_driver_rate=SNV_driver_rate, log_path=log_path)
        setup_logger(file_path=params.log_path)

        if not self._init_clone.library.is_driver_SNV_model:
            params.SNV_driver_rate = 0

        if self.gen == 0:
            logging.info(f'STARTING AGENT GROWTH MODEL')
            logging.info('Site\tNumClones\tPopSize\tMeanFit\tNumBirths\tNumDeaths\tNumMigrations')
            mean_fits = totbirths = totdeaths = totmigrations = [0 for i in range(self.anatomy.nsites)]
            self._log_current_gen(mean_fits, totbirths, totdeaths, totmigrations)
        else:
            logging.info(f'\nCONTINUE')

        TERMINATE = self._terminate(params.t_max, params.min_detectable)
        while not TERMINATE:
            self.gen += 1
            Cell.gen = self.gen
            if self._terminate(params.t_max, params.min_detectable):
                self.gen -= 1
                Cell.gen = self.gen
                TERMINATE = True
                break
            mean_fits, totbirths, totdeaths = self._cycle_birthdeath(params.focal_driver_rate, params.arm_rate, params.chromosomal_rate, params.WGD_rate, params.focal_gain_rate, params.chrom_dup_rate, params.length_mean, params.mag_mean, params.SNV_driver_rate, params.lifespan_mean)
            totmigrations = self._cycle_migrations(params.t_max, pattern=pattern)
            self._log_current_gen(mean_fits, totbirths, totdeaths, totmigrations)
    
    def cycle_gen(
        self, 
        params: Optional[Parameters] = None,
        t_max: Optional[int] = None, 
        lifespan_mean: Optional[Union[int, float]] = None,
        min_detectable: Optional[int] = None, 
        focal_driver_rate: Optional[float] = None, 
        arm_rate: Optional[float] = None, 
        chromosomal_rate: Optional[float] = None, 
        WGD_rate: Optional[float] = None, 
        focal_gain_rate: Optional[float] = None, 
        chrom_dup_rate: Optional[float] = None, 
        length_mean: Optional[Union[int, float]] = None, 
        mag_mean: Optional[Union[int, float]] = None, 
        SNV_driver_rate: Optional[float] = None,
        log_path: Optional[str] = None,
        pattern: str = ''
    ):
        """Identical to the :code:`simulate_agents` method, but runs for a single generation only

        """
        params = fill_params(params, t_max=t_max, lifespan_mean=lifespan_mean, min_detectable=min_detectable, focal_driver_rate=focal_driver_rate, arm_rate=arm_rate, chromosomal_rate=chromosomal_rate, WGD_rate=WGD_rate, focal_gain_rate=focal_gain_rate, chrom_dup_rate=chrom_dup_rate, length_mean=length_mean, mag_mean=mag_mean, SNV_driver_rate=SNV_driver_rate, log_path=log_path)
        setup_logger(file_path=params.log_path)

        if self.gen == 0:
            logging.info(f'STARTING AGENT GROWTH MODEL')
            logging.info('Site\tNumClones\tPopSize\tMeanFit\tNumBirths\tNumDeaths\tNumMigrations')
            mean_fits = totbirths = totdeaths = totmigrations = [0 for i in range(self.anatomy.nsites)]
            self._log_current_gen(mean_fits, totbirths, totdeaths, totmigrations)
        else:
            logging.info(f'\nCONTINUE')
        
        if not self._terminate(params.t_max, params.min_detectable):
            self.gen += 1
            Cell.gen = self.gen
            mean_fits, totbirths, totdeaths = self._cycle_birthdeath(params.focal_driver_rate, params.arm_rate, params.chromosomal_rate, params.WGD_rate, params.focal_gain_rate, params.chrom_dup_rate, params.length_mean, params.mag_mean, params.SNV_driver_rate, params.lifespan_mean)
            totmigrations = self._cycle_migrations(params.t_max, pattern=pattern)
            self._log_current_gen(mean_fits, totbirths, totdeaths, totmigrations)

    def sample_cells(
        self, 
        params: Optional[Parameters] = None,
        ncells_prim: Optional[int] = None,
        ncells_meta: Optional[int] = None, 
    ):
        """Samples cells from each anatomical site. The probability that a cell is sampled from a clone is proportional to that clones population size. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Parameters, optional): 
            ncells_prim (int, optional): 
            ncells_meta (int, optional): 
        """
        params = fill_params(params, ncells_prim=ncells_prim, ncells_meta=ncells_meta)
        obs = [clone for clone in self._observed.keys()]
        for clone in obs:
            del self._observed[clone]
        
        for s in range(self.anatomy.nsites):
            if s == 0:
                ncells = params.ncells_prim
            else:
                ncells = params.ncells_meta
            clones = self.clones[s]
            pops = np.array([clone.popsize for clone in clones])
            tot_cells = sum(pops)
            if ncells > 0 and tot_cells >= ncells:
                idxs = np.random.choice(range(len(clones)), ncells, p=pops/tot_cells)
                for i in idxs:
                    self._observed[clones[i]] += 1
        self._clone_tree = create_clone_tree(self._observed)

    # Master function for creating resolved tree from sampled cells
    def simulate_singlecell_lineage(
        self, 
        params: Optional[Parameters] = None,
        out_dir: Optional[str] = None, 
        focal_pass_rate: Optional[float] = None, 
        SNV_pass_rate: Optional[float] = None, 
        focal_gain_rate: Optional[float] = None, 
        length_mean: Optional[Union[int, float]] = None, 
        mag_mean: Optional[Union[int, float]] = None, 
        bin_size: Optional[int] = None,
        ncells_normal: Optional[int] = None,
        ref: Optional[str] = None,
        alt_ref: Optional[str] = None
    ):
        """Construct a single-cell lineage from the sampled cells under a coalescent, constrained to obey the underlying clonal evolutionary history, and add passenger mutations. Requires running the :code:`sample_cells` class method beforehand. Automatically saves the ground truth tree, migration graph, mutation history, CNA profiles, and SNV profiles. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Params, optional): 
            out_dir (str, optional): 
            focal_pass_rate (float, optional): 
            SNV_pass_rate (float, optional): 
            focal_gain_rate (float, optional): 
            length_mean (int, float, optional): 
            mag_mean (int, float, optional): 
            bin_size (int, optional): 
            ncells_normal (int, optional): 
            ref (str, optional): 
            alt_ref (str, optional): 

        Returns:
            Tree: The single-cell lineage tree.
        """
        params = fill_params(params, out_dir=out_dir, focal_pass_rate=focal_pass_rate, SNV_pass_rate=SNV_pass_rate, focal_gain_rate=focal_gain_rate, length_mean=length_mean, mag_mean=mag_mean, bin_size=bin_size, ncells_normal=ncells_normal, ref=ref, alt_ref=alt_ref)

        if len(self._observed) == 0:
            raise ValueError("Must run GrowthSimulator.sample_cells beforehand.")

        tree = create_singlecell_tree(self._clone_tree, self._observed, self.gen)
        resolve_multifurcation(tree)
        evolve_tree_with_passengers(tree, self._igenome, params.focal_pass_rate, params.SNV_pass_rate, params.focal_gain_rate, params.length_mean, params.mag_mean)
        if params.ncells_normal >= 1:
            add_dummy_diploids(tree, self.gen, params.ncells_normal)
        tree.save(os.path.join(params.out_dir, 'cell_tree_full.nwk'))
        collapse_tree(tree)
        tree.save(os.path.join(params.out_dir, 'cell_tree.nwk'))
        if self.anatomy.nsites > 1:
            edges = generate_migration_graph(tree, self.anatomy.site_ids)
            save_migration_graph(params.out_dir, edges)
        save_mutations_from_tree(tree, params.out_dir, params.bin_size)
        save_CNPs(tree, params.out_dir, params.bin_size)
        save_SNVs(tree, params.out_dir)
        return tree

    def simulate_clonal_lineage(
        self, 
        params: Optional[Parameters] = None,
        out_dir: Optional[str] = None, 
        min_mut_fraction: Optional[float] = None, 
        focal_pass_rate: Optional[float] = None, 
        SNV_pass_rate: Optional[float] = None, 
        focal_gain_rate: Optional[float] = None, 
        length_mean: Optional[Union[int, float]] = None, 
        mag_mean: Optional[Union[int, float]] = None, 
        bin_size: Optional[int] = None,
        ncells_normal: Optional[int] = None,
        ref: Optional[str] = None,
        alt_ref: Optional[str] = None
    ):
        """Construct a clonal lineage from the sampled cells and adds passenger mutations. Optionally prune the observed clones based on a minimum genotype frequency threshold. Requires running the :code:`sample_cells` class method beforehand. Automatically saves the ground truth tree, migration graph, mutation history, CNA profiles, and SNV profiles. See :ref:`Parameters <parameters>` for an explanation of the parameters.

        Args:
            params (Params, optional): 
            out_dir (str, optional): 
            min_mut_fraction (float, optional): 
            focal_pass_rate (float, optional): 
            SNV_pass_rate (float, optional): 
            focal_gain_rate (float, optional): 
            length_mean (int, float, optional): 
            mag_mean (int, float, optional): 
            bin_size (int, optional): 
            ncells_normal (int, optional): 
            ref (str, optional): 
            alt_ref (str, optional): 

        Returns:
            Tree: The clonal lineage tree.
        """
        params = fill_params(params, out_dir=out_dir, min_mut_fraction=min_mut_fraction, focal_pass_rate=focal_pass_rate, SNV_pass_rate=SNV_pass_rate, focal_gain_rate=focal_gain_rate, length_mean=length_mean, mag_mean=mag_mean, bin_size=bin_size, ncells_normal=ncells_normal, ref=ref, alt_ref=alt_ref)

        if len(self._observed) == 0:
            raise ValueError("Must run GrowthSimulator.sample_cells beforehand.")

        tree, obs = create_converted_tree(self._clone_tree, self._observed)
        sample_counts = merge_and_resolve_clone_tree(tree, obs, self.anatomy.nsites, self.gen, params.min_mut_fraction)
        evolve_tree_with_passengers(tree, self._igenome, params.focal_pass_rate, params.SNV_pass_rate, params.focal_gain_rate, params.length_mean, params.mag_mean)
        if params.ncells_normal >= 1:
            add_dummy_diploids(tree, self.gen, 1)
            dip = tree.find('diploid')
            assert dip != None, "Cannot find diploid node..."
            dip.info['count'] = params.ncells_normal
        tree.save(os.path.join(params.out_dir, 'clone_tree_full.nwk'))
        collapse_tree(tree)
        tree.save(os.path.join(params.out_dir, 'clone_tree.nwk'))
        for cell,count in sample_counts.items():
            cell.info['count'] = count
        if self.anatomy.nsites > 1:
            edges = generate_migration_graph(tree, self.anatomy.site_ids)
            save_migration_graph(params.out_dir, edges)
        save_mutations_from_tree(tree, params.out_dir, params.bin_size)
        save_CNPs(tree, params.out_dir, params.bin_size)
        save_SNVs(tree, params.out_dir, ref=params.ref, alt_ref=params.alt_ref)
        site_CN_averages(tree, params.out_dir, params.bin_size)
        return tree
    
    def save_checkpoint(self, outfile):
        """Saves the growth simulator in its current state to a file.

        Args:
            outfile (str): The location to save the growth simulator.
        """
        with open(outfile, 'wb') as f:
            pickle.dump(self, f)
    

