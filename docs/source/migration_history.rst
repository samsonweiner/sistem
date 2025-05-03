Studying distribution of migration patterns
===========================================

In this tutorial, we will show how to create a set of multi-site simulations, and then identify the relative frequency of different migration patterns for the set of parameters used. Migration patterns are determined by the structure of the underlying migration graph :math:`G` based on two criteria. First, :math:`G` is monoclonal (m) if the multiplicity of each multi-edge is 1, and polyclonal (p) otherwise. Second, :math:`G` single-source seeding (S) if :math:`G` is a multi-tree, multi-source seeding (M) if :math:`G` is a multi-DAG, and reseeding (R) otherwise (contains cycles). Thus, any migration graph will either be mS, pS, mM, pM, mR, or pR. Within the mS (resp. pS) patterns, we further subclassify migration graphs in which each site was seeded solely by the primary site. These are designated by mPS (resp. pPS), whereas mS (resp. pS) refer to single-source seeding patterns which contain atleast one metastasis-to-metastasis seedings.

We will be using the Chromosome-Arm selection model with random coefficients and the Static migration model. We will also create site-specific selection libraries which reflect some underlying organotropism priors in the form of a pairwise distance matrix, which will also be randomly generated. The simulations in this example will be for 5 sites, and the per-cell per-generation migration probability between the two farthest sites will be :math:`1e-8`, with the other pairwise migration probabilities scaled according to their distance. The maximum number of coefficients which will differ between any pairs of sites will be :math:`\approx 0.2*44 = 8.8`, as there are 44 coefficients in the Chromosome-Arm model. The agent-based growth phase will terminate after each site reaches a population size of :math:`5e5` cells out of a max carrying capacity of :math:`5e6`, at which point we will sample :math:`10,000` cells from each site and construct a clone tree from the set of observed clones with genotypes occurring at :math:`>=5\%` frequency in the sampled cells of their respective site.

For each simulation, the migration graph will be saved in the respective output directory. To extract the pattern from the graph, we will use some helper functions which you can access `here <https://github.com/samsonweiner/sistem/blob/main/scripts/get_migration_pattern.py>`_.

.. code-block:: python

    import os
    from collections import Counter
    from sistem import Parameters, GrowthSimulator
    from sistem.selection import RandomArmLibrary
    from sistem.anatomy import StaticAnatomy

    # will be using helper functions from get_migration_pattern.py
    from get_migration_pattern import read_graph, get_pattern

    # Set to the desired sample size and where all the simulations will be stored
    num_sims = 1
    sims_dir = 'sim3_out'
    if not os.path.isdir(sims_dir):
        os.makedirs(sims_dir)

    # Create the simulations and save in sims_dir/1, sims_dir/2, ...
    for i in range(num_sims):
        cur_sim_dir = os.path.join(sims_dir, f'{i}')

        # Initialize parameters
        params = Parameters(
            out_dir = cur_sim_dir,
            nsites=5, 
            epsilon=1e-8, 
            min_detectable=5e5, 
            capacities=5e6,
            focal_driver_rate=5e-4, 
            alter_prop=0.2, 
            ncells_prim=10000, 
            ncells_meta=10000, 
            min_mut_fraction=0.05
        )

        # Create a Library instance and initialize with random coefficients
        library = RandomArmLibrary(params=params)
        library.initialize(params=params)

        # Create the Anatomy instance, initialize the distances, and create site-specific libraries.
        anatomy = StaticAnatomy(library, params=params)
        anatomy.initialize_distances()
        anatomy.create_random_metastatic_libraries(method='distance', params=params)

        # Now create the Growth Simulator and start the simulation. This should take around 30 minutes.
        gs = GrowthSimulator(anatomy)
        gs.simulate_agents(params=params)

        # Sample cells and then build the clonal lineage and migration graph. This should only take a few minutes.
        gs.sample_cells(params=params)
        tree = gs.simulate_clonal_lineage(params=params)

        # Now for each simulation we will use some helper functions to identify the pattern.
        patterns = []
        for i in range(num_sims):
            graph_path = os.path.join(sims_dir, f'{i}/migration_graph.tsv')
            g = read_graph(graph_path)
            p = get_pattern(g)
            patterns.append(p)

        # And we can get the relative counts of each pattern as follows:
        c = Counter(patterns)
        for p,count in c.items():
            print(p, count)