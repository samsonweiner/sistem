.. _quickstart:

Quickstart
==========

This page provides a demonstration of a simulation experiment using the main features of SISTEM, with links to more detailed documentation. The main workflow of SISTEM can be broken down into the following key steps:

1. :ref:`Initialize parameters <parameters-label>`
2. :ref:`Initialize selection library <selection-label>`
3. :ref:`Initialize migration model and metastatic landscapes <migration-label>`
4. :ref:`Run agent-based growth simulator <agentgrowth-label>`
5. :ref:`Sample cells and simulate lineage <lineage-label>`
6. :ref:`Generate data <data-label>`

.. _parameters-label:

Initialize parameters
---------------------

SISTEM takes as input numerous user-defined :ref:`parameters <parameters>` in order to generate a wide range of possible mutations, evolutionary histories, migration patterns, and more. When running a simulation, you can pass these parameter values directly on the fly at each step. However, it is advisable that you create an initial set of parameters beforehand. This has the benefit of checking parameter values and can be conveniently passed at each step instead of directly inputting individual parameters.

In this experiment, we are simulating the growth of a human cancer which metastasizes to 2 anatomical sites until detection when each site reaches a population of 1 million cells. For now, don't worry what each parameters mean - these will be explained in the following sections.

.. code-block:: python

    from sistem import Parameters

    params = Parameters(
        out_dir='example_out', 
        nsites=3, 
        epsilon=1e-9, 
        min_detectable=5e5, 
        capacities=1e7, 
        focal_driver_rate=5e-4, 
        SNV_pass_rate=0.01, 
        alter_prop=0.3, 
        ncells_prim=10000, 
        ncells_meta=5000, 
        ncells_normal=2000, 
        min_mut_fraction=0.05,
        coverage=100
    )

The remainder of parameters will be kept to their defaults.

.. _selection-label:

Initialize selection library
----------------------------

SISTEM utilizes an agent-based model to simulate tumor growth and metastasis over discrete generations. The selection library is used to define the relative fitness of a cell based on its genome, which determines the probability of that cell replicating vs dying off and (optionally) successfuly migrating. A complete list of built-in selection models can be found :ref:`here <selectionmodels>`.

In this experiment, we will be using a driver region selection model with random selection coefficients. When initializing a selection library, you can pass some or all of its parameters directly (if unspecified, will revert to defaults) or you can initialize it entirely from a Parameters instance.

.. code-block:: python

    from sistem.selection import RandomRegionLibrary
    from sistem.genome import hg38_chrom_lengths_from_cytoband

    # Create a selection library by directly passing each parameter. Here 
    # we are using chromosome lengths and arm ratios derived from the hg38 human
    # reference genome.
    chrom_lens, arm_ratios = hg38_chrom_lengths_from_cytoband()
    library = RandomRegionLibrary(
        chrom_lens=chrom_lens,
        arm_ratios=arm_ratios,
        region_len=5e6,
        max_ploidy=8,
    )
    # Initialize the selection coefficients
    library.initialize(CN_coeff=0.1, OG_r=0.05, TSG_r=0.05)

    # OR

    # You can parameterize the selection library with just a Parameters instance.
    # By default, the 22 chromosome lengths and arm ratios from hg38 are used.
    library = RandomRegionLibrary(params=params)
    library.initialize(params=params)

Here, :code:`chrom_lens` is a dictionary where keys are chromosome names and values are the length of the chromosome in number of base pairs. You can pass genomes of any size this way, or initialize :code:`chrom_lens` from a reference genome in fasta format. :code:`arm_ratios` is either a dictionary where keys are chromosome names and values are the proportion of the short arm, or a float ratio used for all chromosomes. :code:`region_len` is the size at which the chromosome sequences are broken up into. Lower values increase the resolution of the genome but decrease performance. :code:`max_ploidy` is parameters used in viability checkpoints - if a cell exceeds a ploidy of 8, it immediately dies. :code:`CN_coeff` is the max possible magnitude of selection coefficient for driver genes. :code:`OG_r` and :code:`TSG_r` describe the ratio of regions which are oncogenes and tumor suppressor genes, respectively.

Let's examine the driver regions and their selection coeffients:

.. code-block:: python

    for chrname in params.chrom_lens:
        for region_id in library.drivers[chrname]:
            print(chrname, region_id, library.delta[chrname][region_id])

.. code-block:: console

    chr1 3 -0.07050650925424783
    chr1 23 0.056907147220048596
    chr1 41 -0.01553590386497552
    chr1 48 0.07579778217001987
    chr10 11 -0.058779514551353024
    ...

.. _migration-label:

Initialize migration model and metastatic landscapes
----------------------------------------------------

In this step, we initialize the migration model, which allows us to define organotropism priors and/or genotype-driven migrations, and optionally define site-specific selection libraries. In this experiment, we will use a static model and alter the selection libraries of metastatic sites to reflect these distances. You can find the built-in migration models :ref:`here <migrationmodels>`.

.. code-block:: python

    from sistem.anatomy import StaticAnatomy

    # Pass parameters directly.
    anatomy = StaticAnatomy(
        libraries=library, 
        nsites=3,
        epsilon=1e-8,
        growth_rate=0.005,
        capacities=1e7
    )
    anatomy.initialize_distances(method='random') #method is random or precomputed
    anatomy.create_random_metastatic_libraries(method='distance', alter_prop=0.3, CN_coeff=0.1) #method is random or distance

    # Or pass a Parameters instance.
    anatomy = StaticAnatomy(libraries=library, params=params)
    anatomy.initialize_distances(method='random') #method is random or precomputed
    anatomy.create_random_metastatic_libraries(method='distance', params=params) #method is random or distance

Here, *nsites* is the number of anatomical sites, :code:`epsilon` is the baseline per-generation migration probability, :code:`growth_rate` is the the logistic growth rate, and *capacities* is the carrying capacity of each site (can be an int or a list of ints, one for each site). The :code:`alter_prop` parameter is used for creating metastatic libraries with the 'distance' method and coincides with the ratio of driver genes with a different selection coefficient for the farthest metastatic site from the primary site.

We can view the pairwise site distances (acting as organotropism priors) as a triangular matrix (distances are symmetric).

.. code-block:: python

    import numpy as np

    dist_matrix = np.zeros((anatomy.nsites, anatomy.nsites))
    dist_matrix[np.triu_indices(anatomy.nsites)] = anatomy.dists[:, 0]
    print(dist_matrix)

.. code-block:: console

    [[0.         1.         0.50938716]
     [0.         0.         0.54111119]
     [0.         0.         0.        ]]

.. _agentgrowth-label:

Run agent-based growth simulator
--------------------------------

Now we are ready to run the main stage of the simulator. Starting from an initial clone in the primary site, growth will occur over discrete generations until a either a minimum population size threshold is achieved in each site, all cells die off, or the max number of generations :code:`t_max` is reached. We first create an instance of the *GrowthSimulator* class, then run the *simulate_agents* method. Among other parameters, this function takes a driver focal (segmental) CNA rate :code:`focal_driver_rate` and a driver SNV rate :code:`SNV_driver_rate`, the later of which will be set to 0 because SNVs have no impact on fitness under the driver region selection model. To improve performance, SISTEM only simulates driver mutations (those that alter fitness) during the agent growth stage, with passenger mutations added in the following stage. More information on the growth simulators :ref:`here <growthsim>`.

.. code-block:: python

    import os
    from sistem import GrowthSimulator

    # Initialize GrowthSimulator instance with the anatomy instance.
    gs = GrowthSimulator(anatomy)

    # Run the simulate_agents method
    gs.simulate_agents(
        t_max=5000, 
        min_detectable=5e5, 
        focal_driver_rate=5e-4,
    )

    # Reminder: you can also specify parameters directly with:
    # gs.simulate_agents(params=params)

    # For some experiments, you may want to reuse the simulated populations/history,
    # so its a good idea to create a checkpoint.
    gs.save_checkpoint(os.path.join(params.out_dir, 'gs.pkl'))

We can investigate the clonal populations after completion (or during) during the growth phase. For example, for each site let's print out the number of clones, number of cells, and the fitness of the least and most fit clones after termination.

.. code-block:: python

    for s,clones in gs.clones.items():
        fits = [clone.fitness for clone in clones]
        print(f"Site {gs.anatomy.site_ids[s]}")
        print(f"Num clones: {len(clones)}, Num cells: {gs.site_counts[s]}, Least fit: {min(fits):.2f}, Most fit: {max(fits):.2f}")
        print()

.. code-block:: console
    
    Site P
    Num clones: 4995, Num cells: 9407357, Least fit: 172.64, Most fit: 294.98

    Site A
    Num clones: 230, Num cells: 501285, Least fit: 184.83, Most fit: 284.86

    Site B
    Num clones: 709, Num cells: 1272163, Least fit: 152.55, Most fit: 257.05

.. _lineage-label:

Sample cells and simulate lineage
---------------------------------

Once the agent growth phase completes, the next step involves sampling cells from the wider populations in each site, constructing a lineage (either at the clonal or single-cell level), and adding passenger mutations. In this experiment, we will be generating a clone tree and migration graph from 10,000 cells sampled from the primary site and 5,000 cells sampled from the metastatic sites.

.. code-block:: python

    from sistem import load_gs

    # Load the saved GrowthSimulator instance, if applicable
    gs = load_gs(os.path.join(params.out_dir, 'gs.pkl'))

    # Begin by sampling cell
    gs.sample_cells(
        ncells_prim=10000
        ncells_meta=5000
    )
    # Or run gs.sample_cells(params=params)

    # Call the lineage simulation method
    tree = gs.simulate_clonal_lineage(
        out_dir='example_out', 
        SNV_pass_rate=0.01, 
        ncells_normal=2000,
        min_mut_fraction=0.05
    )
    # Or run gs.simulate_clonal_lineage(params=params)

Here, :code:`out_dir`is the output directly of the experiment, :code:`SNV_pass_rate` is the passenger SNV rate, and *ncells_normal* is the number of additional normal cells sampled from the primary site (a relative number will be sampled from each metastatic site). The last parameter, :code:`min_mut_fraction`, is a way to merge or filter out low-frequency clones from the sample (described more in :ref:`this doc <growthsim>`). 

When running *gs.simulate_clonal_lineage* (or *gs.simulate_singlecell_lineage* if generating a single-cell tree), a number of output data modalities are generated automatically. This includes the ground truth clone tree (or single-cell tree) saved in newick format, the ground truth migration graph, CNA profiles, SNV profiles, and a mutation event log. 


.. _data-label:

Generate count data and sequencing reads
----------------------------------------

Generating read count data or synthetic sequencing reads requires an additional step after the lineage simulation. SISTEM currently can only generate raw single-cell DNA-sequencing reads, so in this experiment we will be generating read count data only. This consists of a total and alternate read count for each SNV occurring in any of the observed clones, with read counts computed separately for each anatomical site. 

.. code-block:: python

    from sistem.data import gen_readcounts_bulk

    # Call the gen_readcounts_bulk function with the tree outputted from the previous
    # step and any non-default parameters.
    gen_readcounts_bulk(
        tree, 
        out_dir='example_out',
        coverage=100
    )
    # Or run gen_readcounts_bulk(tree, params=params)