MPI support
===========

Nonpareil supports MPI (Message Passing Interface) since v2.2. This code is relatively stable, but
it's not as widely tested as the regular Nonpareil.

Requirements
------------

You will first need OpenMPI_ in your computer. There are other MPI implementations, but Nonpareil only supports OpenMPI (by now). Once
you have it, you should have at least the C++ compiler (typically ``mpic++``) and the interactive executable (typically ``mpirun``). If
you have the compiler in a non-standard location (for example, to coexist with mpich), change the value of ``mpicpp`` in the ``globals.mk``
file. Once you are ready, simply run::

   cd nonpareil # or wherever you have the nonpareil folder
   make nonpareil-mpi

That's it. Now you should have the ``nonpareil-mpi`` binary, that you can place in a location listed in your ``$PATH`` if you want.

Running Nonpareil MPI
---------------------

1. Get your machines ready. If you are familiar with MPI skip directly to #3. If you have your own infrastructure, just make sure they
   are MPI-capable (network, permissions, software, etc.). If you are using a cluster, just request as many machines as you need (see
   the resources section below). For example, to request 10 machines with 16 CPUs each in PBS, use ``-l nodes=10:ppn=16``.

2. Obtain the machine names. Just prepare a raw text file with the list of machines you want to use. If you are using PBS, you can do
   this by running::
       
       cat $PBS_NODEFILE | awk 'NR%16==0' > hosts.txt # Change the '16' by the number of CPUs you are using (the value of ppn).

3. Run Nonpareil MPI. All you need is to call ``nonpareil-mpi`` with ``mpirun``. For example, if you want to use 10 machines, with 16
   CPUs each, and the list of machines is in ``hosts.txt``, then run::
       
       mpirun -np 10 -machinefile hosts.txt nonpareil-mpi -t 16 -s path/to/your/sequences.fasta -b output ...

   Note that the options of ``nonpareil-mpi`` are the exact same as for ``nonpareil``. Just remember that the value of ``-t`` is the
   number of threads *per machine*, not the total number of CPUs.

Resources
---------

If you are interested on MPI, I'm assuming you have big files, so you may be also concerned about resources allocation.

How much memory you will need?
   In the `Nonpareil paper`_ (Suppl. Fig. 6) you can see the linear relationship between maximum required RAM and the size of the
   dataset. The function is approximately ``RAM = Size + 2``, where ``RAM`` and ``Size`` are both in Gb. You can use less RAM than
   that, and Nonpareil will adapt, but it'll take longer running times. This value is the "maximum required", which means that if you
   assign more RAM than that, it won't make any difference. Now, that value is the total RAM required. That means that if you use the 
   MPI implementation, you can divide ``Size`` by the number of computers you are using, and then apply the function above. For example,
   if you have a 50Gb dataset, you will need (maximum) 52Gb (50 + 2) of RAM for the standard implementation of Nonpareil. However, if
   you use the MPI version with, for example, 10 machines, you'll need (maximum) 7Gb (50/10 + 2) on each machine.

How many machines you will need?
   I don't have a large benchmarking yet for the MPI version, but at the end it really depends on your resources. If you have more machines,
   it will run faster (unless you have a very small dataset) and it will require less memory (as discussed above).

Should I use more machines or more threads?
   Again, it depends on your resources. Multi-threading is (in general) more efficient, because it doesn't have the overhead of network
   communication. That means that you should favor more CPUs over more machines. However, there are some aspects to take into account. One,
   as discussed above, is the RAM. More machines = less RAM per machine, while more threads have little impact on RAM usage (actually,
   more threads = slighly more RAM). Another catch is the resources availability. It is possible that you have tens of machines for your
   exclusive use, but most likely you are actually sharing resources through a cluster architecture. If you ask for 64 processors per node
   (assuming you have 64-core machines) you will probably have to wait in queue for quite some time. If you ask for 4 machines, and 64
   processors per node, you will likely be waiting in queue for hours or days. However, the same number of threads (256) can be gathered
   by asking for 16 machines, and 16 processors per node. If you do that, you will give the scheduler more flexibility (note that the nodes=4
   ppn=64 is a special case of nodes=16 and ppn=16) hence reducing your queue time. You may be asking: can I simply ask for nodes=256 and ppn=1?
   Well... you can, but as I said multi-threading is more efficient than multi-nodes, so don't go to the extremes. Also, Nonpareil has three
   expensive steps:
   
   1. Reading the fasta, which is strictly linear: only one thread is used in only one machine. This process is linear in time with the size
      of the input file.
   
   2. Comparing reads, which is threaded and multi-node. This is by far the most expensive step, and it is distributed across machines and
      across CPUs on each machine. This process is linear in time with the size of the input file.

   3. Subsampling, which is threaded but not multi-node. This step is not too expensive, and it's nearly constant time. With default parameters,
      it takes about 3 minutes of CPU time, but it grows if you reduce ``-i``. The time on this step is reduced by more threads (``-t``), but
      not by more machines.

How can I evaluate the performance in pilot runs?
   I must say: I rarely do pilot runs. However, I'm often interested on performance for future runs (for example, for other projects). There are
   two sources of information that can be handy. One, is the OS itself (or the PBS output file, if you have a good Epiloge configured). For example,
   to measure the total RAM used, the total walltime, real time, user time, etc. Another source is the .npl file, which contains a log of the
   Nonpareil run (assuming you used the ``-b`` option). The number in squared brackets is the CPU time in minutes. Note that the CPU time here is
   only for the "master" machine. That means: the number of CPU minutes added for all the threads in the main machine. Another useful piece of
   information is the number of "blocks" used. Ideally, you should have one block per machine; if you have more it means that the RAM assigned
   (``-R``) was insufficient. You can find it right below the "Designing the blocks scheme..." line. In the ideal scenario (enough RAM), you should
   have one Qry block, and as many Sbj blocks as machines (one, if you are not using the MPI implementation). If you have more than that, you could
   attain shorter running times by increasing the RAM (``-R``).


.. _OpenMPI: http://www.open-mpi.org/
.. _Nonpareil paper: http://bioinformatics.oxfordjournals.org/content/early/2013/11/05/bioinformatics.btt584.abstract

