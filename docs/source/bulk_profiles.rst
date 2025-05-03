Generating bulk profiles and clone tree
=======================================

This example showcases how to simulate a clonal lineage and directly construct Bulk-seq mutation profiles. The commands used are very similar to the tutorial described in the :ref:`quickstart <quickstart>` section. However, this time we will be using the Chromosome-Arm selection model with random coefficients and the Genotype migration model. Running these commands will automatically generate the clonal lineage tree, migration graph, and bulk mutation profiles. In particular, it will generate clone-level CNA and SNV profiles and site-wise copy number averages.

The simulations in this example will be for 3 sites and a baseline per-cell per-generation migration probability of :math:`1e-8`. Under the genotype model, the probability that cell :math:`c` in site :math:`a` migrates to site :math:`b` is at minimum :math:`1e-8` but increases as :math:`c` becomes more fit with respect to the selection library of :math:`b`. The other parameters will be the same as in the quickstart tutorial.

.. code-block:: python

    import os
    from sistem import Parameters, GrowthSimulator
    from sistem.selection import RandomArmLibrary
    from sistem.anatomy import GenotypeAnatomy

    params = Parameters(
        out_dir='sim2_out', 
        nsites=3, 
        epsilon=1e-8, 
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

    # Create a Library instance and initialize with random coefficients
    library = RandomArmLibrary(params=params)
    library.initialize(params=params)

    # Create the Anatomy instance and initialize the distances.
    anatomy = GenotypeAnatomy(library, params=params)
    anatomy.initialize_distances()

    # Now create the GrowthSimulator.
    gs = GrowthSimulator(anatomy)

    # Start the simulation. This may take ~20 minutes. You can monitor the sim.log file in real time.
    gs.simulate_agents(params=params)
    gs.save_checkpoint(os.path.join(params.out_dir, 'gs.pkl'))
    gs.sample_cells(params=params)
    tree = gs.simulate_clonal_lineage(params=params)

