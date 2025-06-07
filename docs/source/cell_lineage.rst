.. _scCNAtree:

Generating single-cell CNA profiles and cell tree
=================================================

In this example, we will demonstrate how to use SISTEM to generate CNA profiles from observed cells with a known ground truth cell lineage tree. Here we will grow a single tumor until it reaches 5,000,000 cells with higher slightly higher chromosome-arm and whole-chromosomal CNA rates than the defaults. We will use the region selection model with random coefficients, and because we are simulating only one site, use the SimpleAnatomy class. Because we are only interested in the CNAs, we will disable SNVs by setting the passenger SNV rate to 0 (the region selection model automatically sets the driver SNV rate to 0).

Once the tumor is simulated, we will create 10 replicate trees by independently sampling 1000 cells and adding random passenger mutations. This will lead to different cell lineage trees and CNA profiles from the same simulated tumor. Note that the simulate_singlecell_lineage method will automatically generate the ground truth tree, migration graph, mutation history, CNA profiles, and SNV profiles. This example may take around 25 minutes to execute.

.. code-block:: python

    import os
    from sistem import GrowthSimulator, Parameters, load_gs
    from sistem.selection import RandomRegionLibrary
    from sistem.anatomy import SimpleAnatomy

    # Initialize parameters
    params = Parameters(
        out_dir='sim1_out',
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
        ncells_prim = 1000
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

    # Now for nreplicates, we will create a subdirectory and simulate a unique cell-lineage tree.
    nreplicates = 10
    for i in range(nreplicates):
        gs = load_gs(os.path.join(params.out_dir, 'gs.pkl'))
        gs.sample_cells(params=params)
        gs.simulate_singlecell_lineage(params=params, out_dir=os.path.join(params.out_dir, f"rep{i}"))