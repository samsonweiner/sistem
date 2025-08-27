.. SISTEM documentation master file, created by
   sphinx-quickstart on Tue Apr 29 23:23:23 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SISTEM documentation
====================

**SISTEM** is an open-source Python package for simulating tumor growth, metastasis, and DNA-seq data under genotype-driven selection. It can be used to generate allele-specific copy number aberration (CNA) profiles, single-nucleotide variant (SNV) profiles, bulk-sequencing read counts, single-cell allele-specific read counts, and single-cell DNA sequencing reads, in addition to ground truth clone trees, single-cell lineage trees, and migration graphs. This makes **SISTEM** ideal for benchmarking methods for mutation detection, phylogenetic inference, and migration history inference. Additionally, **SISTEM** provides a variety of customizable selection models and migration models, making it a suitable framework for other applications such as estimating parameters from real datasets.

This manual offers instructions for :ref:`installing <install>` and using **SISTEM** and contains example workflows. You can find an overview of the key steps in the :ref:`quickstart <quickstart>` guide and a complete list of parameters :ref:`here <parameters>`. SISTEM is implemented as an API for high customization, but users can also use the `SISTEM wrapper <https://github.com/samsonweiner/sistem/tree/main/wrapper>`_ for easy end-to-end simulation workflows. If you would like to learn more about the underlying algorithms used in **SISTEM**, please see our paper (yet to be posted).

.. `Bioinformatics paper <https://www.youtube.com>`_.

**SISTEM** can be cited as follows:

| SISTEM: simulation of tumor evolution, metastasis, and DNA-seq data under genotype-driven selection
| Samson Weiner and Mukul S. Bansal
| Under Review

Contents:
=========

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   install
   quickstart

.. toctree::
   :maxdepth: 1
   :caption: API core and features

   selection
   migration
   growth_sim
   seq_count_data
   parameters
   output_formats.rst

.. toctree::
   :maxdepth: 1
   :caption: Tutorials and examples

   cell_lineage
   bulk_profiles
   migration_history
   sequencing_reads