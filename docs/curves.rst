Nonpareil curves
================

The estimation of the :doc:`redundancy` is at the core of Nonpareil, but it's
when those values are transformed into average coverage that they become
comporable across samples, and become useful for project design and sample
evaluation.

To build Nonpareil curves, you need two things. First, the Nonpareil.R file
(you can find it in the ``utils`` folder of Nonpareil). Second, the ``.npo``
file (or ``-o`` value, if you used this option) generated in the estimation of
:doc:`redundancy`.

For the impatient
-----------------

First, load the package. If you don't have it installed yet, you can open R_ and
execute::
   
   install.packages('Nonpareil');
   library(Nonpareil);

If you did `make install` (:doc:`installation`), you can simply open R_ and
execute::

   library(Nonpareil);

And you can get help messages using any of::

   ?Nonpareil.curve
   ?Nonpareil.set
   ?Nonpareil.legend
   ?Nonpareil.predict

Now, you can simply execute::

   Nonpareil.curve('output.npo'); # Change output.npo to the actual redundancy file.

Nonpareil.curve()
-----------------

This function can generate a Nonpareil curve from a ``.npo`` file. See the
documentation of this function inside R_ after loading the Nonpareil package::

   ?Nonpareil.curve
   
Nonpareil.set()
-----------------------

This function can generate a plot with several Nonpareil curves from ``.npo``
files. See the documentation of this function in R_ after loading the Nonpareil
package::

   ?Nonpareil.set

**Example**: I find it very convenient to first prepare a table with the
samples, something like::

    # samples.txt
    File	Name	R	G	B
    SRS063417.1.L50.npo	Posterior fornix	255	200	200
    SRS063287.1.L50.npo	Buccal mucosa	255	120	120
    SRS062540.1.L50.npo	Tongue dorsum	255	3	3
    SRS016335.1.L50.npo	Stool	200	135	76
    SRS015574.1.L50.npo	Supragingival plaque	230	100	120
    SRS019087.1.L50.npo	Anterior nares	220	220	130

Note that this table is tab-delimited, because I find it easier to read, but you
can use anything you like (and is supported by R_). Next, you can simply type
something like this in the R_ console::

    library(Nonpareil);
    samples <- read.table('samples.txt', sep='\t', header=TRUE);
    attach(samples);
    nps <- Nonpareil.set(File, col=rgb(R,G,B,maxColorValue=255), labels=Name,
          plot.opts=list(plot.observed=FALSE));
    detach(samples);
    summary(nps);

To execute examples with real data included in the package, you can execute::

   example(Nonpareil.curve);
   example(Nonpareil.set);

.. _R: http://www.r-project.org/

