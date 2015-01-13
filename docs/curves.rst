Nonpareil curves
================

The estimation of the :doc:`redundancy` is at the core of Nonpareil, but it's when those values are transformed
into average coverage that they become comporable across samples, and become useful for project design and sample
evaluation.

To build Nonpareil curves, you need two things. First, the Nonpareil.R file (you can find it in the ``utils`` folder
of Nonpareil). Second, the ``.npo`` file (or ``-o`` value, if you used this option) generated in the estimation of
:doc:`redundancy`.

For the impatient
-----------------

First, load the package. If you did `make install` (:doc:`installation`), you can simply open R_ and execute::

   library(Nonpareil);

And you can get help messages using any of::

   ?Nonpareil.curve
   ?Nonpareil.curve.batch
   ?Nonpareil.legend
   ?Nonpareil.f
   ?Nonpareil.antif
   ?Nonpareil.coverageFactor

If you didn't install it, you have to load it from the source (although you won't have the embedded documentation)::
   
   source('utils/Nonpareil.R'); # Change utils/Nonpareil.R for the actual path to the utils folder

Now, you can simply execute::

   Nonpareil.curve('output.npo'); # Change output.npo to the actual redundancy file.

Nonpareil.curve()
-----------------

This function can generate a Nonpareil curve from a ``.npo`` file. See the documentation of this function inside R_ after
loading the Nonpareil package::

   ?Nonpareil.curve
   
If you didn't install the Nonpareil package, you can see the documentation from the source::

   tools::Rd2txt(tools::parse_Rd('utils/nonpareil/man/Nonpareil.curve.Rd'))


Nonpareil.legend()
------------------

This function creates a legend for the Nonpareil curve(s) in the (active) plot. It's compatible with single
or multiple calls of `Nonpareil.curve()`_ (using ``new=F`` in all but the first call) and with
`Nonpareil.curve.batch()`_. See the documentation inside R_ after loading the Nonpareil package::

   ?Nonpareil.legend

Or from the source::

   tools::Rd2txt(tools::parse_Rd('utils/nonpareil/man/Nonpareil.legend.Rd'))

Nonpareil.curve.batch()
-----------------------

This function can generate a plot with several Nonpareil curves from ``.npo`` files. See the documentation of this
function in R_ after loading the Nonpareil package::

   ?Nonpareil.curve.batch

Or from the source::

   tools::Rd2txt(tools::parse_Rd('utils/nonpareil/man/Nonpareil.curve.batch.Rd'))

**Example**: I find it very convenient to first prepare a table with the samples, something like::

    # samples.txt
    File	Name	R	G	B
    # HMP
    SRS063417.1.L50.npo	Posterior fornix	256	200	200
    SRS063287.1.L50.npo	Buccal mucosa	256	120	120
    SRS062540.1.L50.npo	Tongue dorsum	256	3	3
    SRS016335.1.L50.npo	Stool	200	135	76
    SRS015574.1.L50.npo	Supragingival plaque	230	100	120
    SRS019087.1.L50.npo	Anterior nares	220	220	130

Note that this table is tab-delimited, because I find it easier to read, but you can use anything you like (and is
supported by R_). Next, you can simply type something like this in the R_ console::

    library(Nonpareil); # Or source('utils/Nonpareil.R');, if you didn't "make install"
    samples <- read.table('samples.txt', sep='\t', h=T);
    attach(samples);
    np <- Nonpareil.curve.batch(File, r=R, g=G, b=B, libnames=Name, modelOnly=TRUE);
    Nonpareil.legend('bottomright');
    detach(samples);


.. _R: http://www.r-project.org/

