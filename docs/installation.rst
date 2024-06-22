Install Nonpareil
====================

Nonpareil can be installed using
`conda <https://bioconda.github.io/recipes/nonpareil/README.html>`_, as a
`Biocontainer <https://quay.io/repository/biocontainers/nonpareil>`_, with
`Galaxy <https://galaxyproject.org/>`_, through `Homebrew <https://brew.sh>`
or via the source code.

Conda installation
------------------

1. Install `Miniconda <https://conda.io/miniconda.html>`_
2. Configure the channels to access `Bioconda <https://bioconda.github.io>`_::

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

3. Install Nonpareil::

    conda install nonpareil

Biocontainer
------------

1. Install `Docker <https://docs.docker.com/engine/installation/>`_
2. Pull the container::

    $ docker pull quay.io/biocontainers/nonpareil

3. Launch the container::

    $ docker run -i -t quay.io/biocontainers/nonpareil /bin/bash

Galaxy
------

You can install Nonpareil on your own Galaxy instance:

1. Go the Galaxy admin space
2. Search on the main `Toolshed <https://toolshed.g2.bx.psu.edu/>`_ for the
   nonpareil repository available under the "Metagenomics" sections
3. Install it

It will automatically install Nonpareil via the conda installation

Homebrew
--------

You can install Nonpareil using `Homebrew <https://brew.sh>` or
`Linuxbrew <http://linuxbrew.sh/>`.

1. Install `Homebrew <https://brew.sh>` if you haven't yet::

    $ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

2. Get the `brewsci/bio <https://brewsci.github.io/homebrew-bio/>` tap if you haven't yet::

    $ brew tap brewsci/bio

3. Install Nonpareil::

    $ brew install nonpareil

Source code installation
------------------------

System requirements
*******************

**Nonpareil binary**: Nonpareil requires a C++ compiler. It has been tested on
64-bit machines with GCC versions >=4.2.1, running Mac OSX and Red Hat Linux.

**Nonpareil MPI**: If you want to compile Nonpareil with MPI support, you will
need OpenMPI_ (v > 1.4.3 tested). Other implementations of MPI could work, but
are yet untested.

**Nonpareil utilities**: Requires R_. No additional libraries are necessary.

Compilation
***********

1. **Get the source**

   Clone the repository from GitHub_::

      git clone git://github.com/lmrodriguezr/nonpareil.git

   If you don't have git_, you can also download the TAR-Ball_ and unpack it
   with::

      tar zxvf nonpareil.tar.gz

2. **Compile**

   Change directory into the newly created folder, and compile Nonpareil::

      cd nonpareil
      make

   If you want to compile Nonpareil MPI (see also :doc:`mpi`), just run::

      make nonpareil-mpi

   In either case, you can specify the C++ compiler to be used setting the
   ``cpp`` or ``mpicpp`` variables, respectively. For example::

      # This compiles nonpareil with /usr/local/bin/g++
      make cpp=/usr/local/bin/g++ nonpareil
      # This compiles nonpareil-mpi with /usr/local/bin/mpic++
      make mpicpp=/usr/local/bin/mpic++ nonpareil-mpi

3. **Install**

   If you want to make Nonpareil available system-wide, just run::

      sudo make install

   If you don't have superuser privileges and/or want to install Nonpareil in a
   location other than ``/usr/local``, simply set the prefix, for example::

      make prefix=$HOME/apps install

   You can also change the location of ``R``, if it's not in the ``$PATH`` or
   you want to use a non-standard installation::

      make prefix=$HOME R=~/bin/R install

   Other variables you can set explicitly for the ``install`` target are
   ``bindir`` (binaries directory) and ``mandir`` (documentation directory).


.. _R: http://www.r-project.org/
.. _git: http://git-scm.com/
.. _GitHub: https://github.com/lmrodriguezr/nonpareil
.. _OpenMPI: http://www.open-mpi.org/
.. _TAR-Ball: https://github.com/lmrodriguezr/nonpareil/tarball/master
