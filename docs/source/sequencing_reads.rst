Generating synthetic scDNA-seq reads
====================================

This tutorial illustrates how to generate raw synthetic single-cell DNA-sequencing reads using SISTEM. Our experimental setup will mostly follow the :ref:`Generating single-cell CNA profiles and cell tree <scCNAtree>` tutorial, but we will only create one replicate and generating the reads requires an extra step. In addition to the scDNA-seq reads, it will also generate ground truth haplotype-specific CNA profiles and haplotype-specific read counts. As this is a single-cell experiment, we will be using a relatively low mean coverage of :math:`0.02`, and a read length of :math:`150`. To speed up the read simulation, we wil leverage multiple processors. The main simulation step may take around 25 minutes, while generating the sequencing reads will take much longer.

.. code-block:: python

    import os
    from sistem import GrowthSimulator, Parameters
    from sistem.selection import RandomRegionLibrary
    from sistem.anatomy import SimpleAnatomy
    from sistem.data import gen_reads

    # First initialize parameters
    params = Parameters(
        out_dir='sim4_out',
        nsites=1,
        capacities=1e7,
        min_detectable=5e6,
        arm_rate=5e-4,
        chromosomal_rate=1e-4,
        SNV_pass_rate=0,
        region_len=1e6,
        CN_coeff=0.1,
        OG_r=0.1,
        TSG_r=0.1,
        ncells_prim = 1000,
        coverage = 0.02,
        read_len = 150,
        num_processors = 4
    )

    # Create a Library instance and initialize with random coefficients
    library = RandomRegionLibrary(params=params)
    library.initialize(params=params)

    # Create the Anatomy instance and then GrowthSimulator. No need to initialize distances using SimpleAnatomy.
    anatomy = SimpleAnatomy(library, params=params)
    gs = GrowthSimulator(anatomy)

    # Simulate tumor growth. This may take ~20 minutes. You can monitor the sim.log file in real time.
    gs.simulate_agents(params=params)
    gs.save_checkpoint(os.path.join(params.out_dir, 'gs.pkl'))

    # Now sample cells and simulate the cell-lineage tree
    gs.sample_cells(params=params)
    tree = gs.simulate_singlecell_lineage(params=params)
    
    # Call the gen_reads function with the first parameter being the Tree object outputted by gs.simulate_singlecell_lineage.
    gen_reads(tree, params=params)

