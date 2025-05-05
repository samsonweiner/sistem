.. _install:

Installation
====================

Recommended installation
------------------------

You can install SISTEM using :code:`pip` as follows:

.. code-block:: bash

        pip install sistem

While this is sufficient for most use cases, if the intended use is to generate synthetic DNA sequencing reads, you must also have the external packages :code:`samtools` and :code:`dwgsim` installed with the binaries added to your :code:`$PATH` variable. These can easily be installed with bioconda:

.. code-block:: bash

    conda install -c bioconda samtools dwgsim

OSX users may experience difficulties installing the *dwgsim* package with bioconda. If this is the case, or other issues arise, you can easily install the binaries manually (see the following section).


..
    SISTEM can be easily installed with :code:`conda` from the `bioconda <https://bioconda.github.io/>`_ channel. If you have not used bioconda before, run the one-time setup command:

    .. code-block:: bash

        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda config --set channel_priority strict

    For best practices, install SISTEM into a new environment. You can do so as follows:

    .. code-block:: bash

        conda create -n SISTEM
        conda activate SISTEM
        conda install -c bioconda sistem

    Make sure to activate the SISTEM environment before every use. 

    If installing with conda on a OSX machine, you may experience difficulties installing the *dwgsim* package. This package is used to generate synthetic sequencing reads, and is only necessary if that is the intended use of SISTEM. If this is the case, you can install a binary for dwgsim manually from `here <https://github.com/nh13/DWGSIM/blob/main/docs/02_Installation.md>`_. For all other applications, you need only install the required python package dependencies, which are listed under the following section.

    You can also install SISTEM using pip:

    .. code-block:: bash

        pip install sistem

    However, if the intended use is to generate synthetic sequencing reads, you must also manually install :code:`samtools` and :code:`dwgsim` binaries (see below).

Manual installation
-------------------
Users can install SISTEM manually by downloading the source code and setting up an environment with all the necessary packages. The remainder of this section provides instructions on how to install SISTEM this way.

You can download the source code by cloning this repository:

.. code-block:: bash

    git clone https://github.com/samsonweiner/sistem.git


Then, cd into the SISTEM repository, make sure setuptools is installed, and install with pip:

.. code-block:: bash

    pip install setuptools
    pip install .

SISTEM depends on the following Python packages:

* `Numpy <https://numpy.org/>`_
* `Scipy <https://scipy.org/>`_
* `scikit-learn <https://scikit-learn.org/stable/>`_
* `msprime <hhttps://tskit.dev/msprime/docs/latest/intro.html>`_
* `Biopython <https://biopython.org/>`_
* `Pyfaidx <https://github.com/mdshw5/pyfaidx>`_

Python packages can be easily installed with a package manager such as :code:`pip` or :code:`conda`. If using :code:`pip`, you can install the packages by running:

.. code-block:: bash

    pip install numpy scipy msprime biopython pyfaidx

If the intended use of SISTEM is to generate synthetic sequencing reads, you must also have the following external packages be installed and configured on the environment:

* `samtools <http://www.htslib.org/download/>`_
* `dwgsim <https://github.com/nh13/DWGSIM>`_ version >=0.1.13

You can also install the packages individually by following the instructions found on the package homepage. If this option is used, the downloaded binaries must be compiled and the directories containing each binary must be added to your :code:`$PATH` variable. For example,

.. code-block:: bash

    export PATH=/path/to/samtools/bin:$PATH

You may also wish to add this line to your *~/.bashrc* or */.bash_profile* configuration file to avoid having to retype this command on login. 