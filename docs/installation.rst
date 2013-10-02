Install Nonpareil
====================

System requirements
-------------------

**Nonpareil binary**: Nonpareil requires a C++ compiler. It has been tested on 64-bit machines with GCC versions ≥4.2.1, running Mac OSX and Red Hat Linux.

**Nonpareil utilities**: Requires R_. No additional libraries are necessary.

Compilation
-----------

1. **Get the source**

   Clone the repository from GitHub_::

      git clone git://github.com/lmrodriguezr/nonpareil.git

   If you don't have git_, you can also download the TAR-Ball_ and unpack it with::

      tar zxvf nonpareil.tar.gz

2. **Compile and install**

   Change directory into the newly created folder, and compile Nonpareil::

      cd nonpareil
      make

   If you want to make Nonpareil available system-wide, just copy the generated binary into a folder listed in the ``$PATH``. For example::

      sudo cp nonpareil /usr/local/bin/

.. _R: http://www.r-project.org/
.. _git: http://git-scm.com/
.. _GitHub: https://github.com/lmrodriguezr/nonpareil
.. _TAR-Ball: https://github.com/lmrodriguezr/nonpareil/tarball/master
