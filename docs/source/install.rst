.. _install:

Installation
====================

Recommended installation
------------------------
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

**Note:** If installing using conda, you may experience difficulties installing a dependency package named *dwgsim*, particularly on OSX machines. You can install SISTEM without dwgsim for Python version >=3.9 using pip:

.. code-block:: bash

    pip install sistem


The dwgsim package is used to generate synthetic sequencing reads, and is only necessary if that is the intended use of SISTEM. For all other applications, installing SISTEM with pip as-is will be sufficient. If users experience difficulties installing with conda but intend to generate sequencing reads, you can first install SISTEM using pip and then install a binary for dwgsim manually from `here <https://github.com/nh13/DWGSIM/blob/main/docs/02_Installation.md>`_. You will also need an executable for `samtools <http://www.htslib.org/download/>`_. Please see additional details on installing these two dependencies in the following section.


Manual installation
-------------------
Users can install SISTEM manually by downloading the source code and setting up an environment with all the necessary packages. The remainder of this section provides instructions on how to install SISTEM this way.

You can download the source code by cloning this repository:

.. code-block:: bash

    git clone https://github.com/samsonweiner/sistem.git

SISTEM depends on the following Python packages, which can be installed with package managers like :code:`pip` or :code:`conda`.

* `Numpy <https://numpy.org/>`_
* `Scipy <https://scipy.org/>`_
* `scikit-learn <https://scikit-learn.org/stable/>`_
* `msprime <hhttps://tskit.dev/msprime/docs/latest/intro.html>`_
* `Biopython <https://biopython.org/>`_
* `Pyfaidx <https://github.com/mdshw5/pyfaidx>`_
* `setuptools <https://pypi.org/project/setuptools/>`_

Once installed, cd into the SISTEM repository and install the project:

.. code-block:: bash

    pip install .


If generating synthetic sequencing reads, SISTEM also requires that the following external packages be installed and configured on the environment.

* `samtools <http://www.htslib.org/download/>`_
* `dwgsim <https://github.com/nh13/DWGSIM>`_ version >=0.1.13

You can also install the packages individually by following the instructions found on the package homepage. If this option is used, the downloaded binaries must be compiled and the directories containing each binary must be added to your :code:`$PATH` variable. For example,

.. code-block:: bash

    export PATH=/path/to/bin:$PATH

You may also wish to add this line to your *~/.bashrc* or */.bash_profile* configuration file to avoid having to retype this command on login. 