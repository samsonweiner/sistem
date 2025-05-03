.. _outputformats:

Output file information
=======================

sim.log
--------
This file contains a record of the clone count, cell count, and fitness of each site at every generation.

cell_tree.nwk
--------------
The single-cell lineage tree in newick format. There is also a cell_tree_full.nwk which is the same tree but includes redundant intermediate common ancestors. Created with GrowthSimulator.simulate_singlecell_lineage.

clone_tree.nwk
---------------
The clonal lineage tree in newick format. There is also a clone_tree_full.nwk which is the same tree but includes redundant intermediate common ancestors. Created with GrowthSimulator.simulate_clonal_lineage.

migration_graph.tsv
-------------------
The migration graph describing the set of migrations present in the cell or clone tree. Each row has two tab-seperated columns and describes a multi-edge in the graph. The first column is of the form 'X-Y' and indicates the presence of a migration from site X to site Y. The second column is an integer describing the multiplicity of the multi-edge, i.e. the number of independent migrations between these sites. Only created for multi-site simulations.

CNA_events.tsv
---------------
A history of all CNA events occurring in the sampled cells and their ancestors. Columns are Cell/Clone, Type, Driver, Chrom, Allele, Arm, Start, Index, Length, Copies, Ref, Region Indices. 'Start' refers to the stating reference region id of the event, whereas 'Index' is the index of that region in the sequence array. The full list mutated reference regions appear in the 'Region Indices' column.

SNV_events.tsv
---------------
A history of all SNV events occurring in the sampled cells and their ancestors. Columns are Cell/Clone, Type, Driver, Chrom, Allele, Region Index, Ref Region Index, Position, BP. Here 'Region Index' is the index of that region in the sequence array, whereas 'Ref Region Index' is the reference region id. Position is the base pair within the region where the SNV has occurred (will be between 0 and *region_len*). BP is either 0, 1, or 2, and corresponds to the index of the alternate base pair from among [A,C,G,T]/X, where X is the reference base pair. For example, if BP = 1 and X = 'C', then BP = [A,G,T][1] = 'G'.

observed_CNPs.tsv
------------------
The copy number profiles (CNPs) of all observed cells or clones in the tree. Columns are Cell/Clone, Chrom, Start, End, Hap CN, Total CN. Here, Hap CN contains the haplotype-specific copy numbers :math:`x_a`, :math:`x_b`, whereas Total CN contains the total :math:`x = x_a + x_b`.

ancestral_CNPs.tsv
-------------------
Same as observed_CNPs.tsv, but contains the CNPs of all ancestral cells.

SNV_profiles.tsv
-----------------
The genotypes of all observed cells or clones in the tree at each SNV position occurring in any cell/clone. Columns are Cell/Clone, Chrom, Pos, GT. 

readcounts.tsv
---------------
The allele-specific readcounts of the observed cells or clones in the tree. Follows a different format if generated for single-cells with sistem.data.gen_readcounts_singlecell vs for clones with sistem.data.gen_readcounts_bulk. If the former, follows the same format as observed_CNPs.tsv but instead of Hap CN and Total CN columns, there are Acount and Bcount columns containing the read count for allele A and B, respectively. If the latter, the columns are Chrom, Pos, P, A, B, ... P is the primary site, and each metastatic site is given its own column A, B, etc. The values under site column X are a pair :math:`r_t`, :math:`r_v`, where :math:`r_t` is the total read count at that position from that site sample and :math:`r_v` is the variant read count.

CN_averages.tsv
----------------
Contains the site-wide copy number averages, computed as the weighted average of the copy number of all observed clones in the tree, with weights equal to number of sampled cells from that clone. Created with GrowthSimulator.simulate_clonal_lineage. Columns are Chrom, Start, End, Site, CN, where CN is the total copy number average.