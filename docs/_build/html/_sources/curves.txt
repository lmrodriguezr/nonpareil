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

You can simply open R_ and execute::

    source('utils/Nonpareil.R');
    Nonpareil.curve('output.npo');

Changing ``utils/Nonpareil.R`` and ``output.npo`` to point to the actual location of the files.

Nonpareil.curve()
-----------------

This function can generate a Nonpareil curve from a ``.npo`` file. The parameters that this function accepts are:

file
   Path to the ``.npo`` file, containing the read redundancy.

overlap=NULL
   Value of the ``-L`` parameter (in Nonpareil, the default is ``50``). If not set, it tries to find the value in the
   ``.npo`` file (supported in Nonpareil ≥2.0), or fails with an error message.

factor=1
   Factor by which the number of reads must be multiplied (x-axis). For example, a factor of 1e-3 will produce curves
   in Kbp (instead of bp). This option should be used when if an "integer overload" is produced, at risk of reducing
   the accuracy of the curves.

plotDispersion=NA
   Indicates if and how the dispersion of replicates is to be represented. If it's ``NA`` no dispersion is plotted.
   Otherwise, it can take the values ``'sd'``: 1 standard deviation around the mean; ``'ci95'``: 95% confidence
   interval; ``'ci90'``: 90% confidence interval; ``'ci50'``: 50% confidence interval; ``'iq'``: Inter-quartile
   range.

xmin=1e3, xmax=10e12
   Range to plot in the x-axis (Sequencing effort).

ymax=1, ymin=1e-6
   Range to plot in the y-axis (Average coverage).

xlab=NULL,ylab=NULL
   Labels of the axes. If ``NULL``, the default values are set Sequencing effort (*units*) for the x-axis, and
   Estimated average coverage for the y-axis.

r=NA,g=NA,b=NA,
   Values for the red, green, and blue components of the rgb representation of the curve's color. If ``NA``, random
   values are set. If the numbers are in the range [0,1], they are assumed to be fractions of 1 (as interpreted by
   rgb). Otherwise, they are assumed to be fractions of 256 (as represented in most color palettes).

new=TRUE
   Indicates if a new canvas must be created. If ``FALSE``, it requires that Nonpareil.curve() was previously called
   (so a proper canvas exists already). This is useful to create a single plot with several samples (see also
   `Nonpareil.legend()`_).

plot=TRUE
   Indicates if the curve is to be plotted. If ``FALSE``, no plots are generated, but the estimations are.

libname=NA
   Name of the library (to be used in later calls of `Nonpareil.legend()`_). If ``NA`` the base of the filename
   (withouth .npo) is used.

modelOnly=FALSE
   If ``TRUE``, prints only the fitted model. It also causes variables related to dispersion and curve to be
   ignored, and it produces an empty circle denoting the actual size of the dataset and the estimated average
   coverage.

plotModel=TRUE
   Indicates whether the fitted model is to be included in the plot. If ``FALSE``, it still fits the model to
   find the estimate of sequencing effort, but the curve of the model is not included in the plot.

curve.lwd=2
   Line width of the curve (the observed curve, not the fitted model).

curve.alpha=0.4
   Transparency of the curve.

model.lwd=1
   Line width of the fitted model line.

model.alpha=1
   Transparency of the fitted model line.

log='x'
   Use logarithmic scale in the X or Y-axis (or both).

data.consistency=TRUE
   Check consistency of the input data before plotting and fitting the model.

useValue='mean'
   Value of the distribution of replicates to be used as the curve. The fitted model and the estimated coverage
   is only meaningful when ``useValue='mean'``, but all of the following values are supported to allow visual
   inspection of the distribution of replicates: The mean (``'mean'``) or the median (``'median'``) of the
   distribution (``'median'``); the upper (``'ub'``) or the lower (``'lb'``) bound of the 95% confidence interval
   (assuming normality); or the first (``'q1'``) or the third (``'q3'``) quartile of the ditribution.

star=95
   Value of average coverage (in percentage) to be used as *almost complete covereage*. Usually this value should
   be set to 95 or 99, but other estimations can be useful. For example, setting this value to 60 predicts the
   minimum sequencing effort at which an assembly with N50 larger than 200bp can be achieved with (100bp-long)
   Illumina reads.

read.length=101
   Length of the reads. This value can be found in the ``.npl`` file, but other estimation (*e.g.*, the average
   length after trimming) can be used.

**Value**: A list with the following keys:

kappa
   Redundancy of the dataset, as calculated by Nonpareil.

C
   Estimated abundance-weighted average coverage.

LRstar
   Estimated sequecing effort (in bp) required to reach ``star`` average coverage (95%, by default).

LR
   Size of the datasets (in bp).

modelR
   Pearson's correlation coeficient between the calculated values and the fitted model.


Nonpareil.legend()
------------------

This function creates a legend for the Nonpareil curve(s) in the (active) plot. It's compatible with single
or multiple calls of `Nonpareil.curve()`_ (using ``new=F`` in all but the first call) and with
`Nonpareil.curve.batch()`_. The parameters that this function accepts are:

x=NULL
   Position in the X-axis. If ``NULL``, it's located at 75% of the maximum value of X. It can also be set to
   any character string supported by xy.coords.

y=.3
   Position in the Y-axis.

...
   Any other parameter accepted by ``legend()`` is supported, except for ``fill`` and ``legend``.

Nonpareil.curve.batch()
-----------------------

This function can generate a plot with several Nonpareil curves from ``.npo`` files. The parameters that this
function accepts are:

files
   Vector of characters with the paths to the ``.npo`` files.

overlap=NULL
   Value of the ``-L`` parameter (in Nonpareil, the default is ``50``). It can be a number (if all the curves were
   generated with the same value) or a vector (in the same order of ``files``). See the ``overlap`` value of
   `Nonpareil.curve()`_.

r=NA,g=NA,b=NA
   Values of the corresponding ``r``, ``g``, and ``b`` in `Nonpareil.curve()`_. It can be ``NA`` (to set random
   colors) or vectors of numbers in the same order of ``files``.

libnames=NA
   Vector of names of the libraries (corresponding to ``libname`` in `Nonpareil.curve()`_). It must be characters,
   not factors.

read.lengths=NA
   A vector of numbers indicating the length of the reads (corresponding to ``read.length`` in `Nonpareil.curve()`_).
   If ``NA``, the default is used.

...
   Any other parameter accepted by `Nonpareil.curve()`_ is supported.

**Value**: A dataframe containing the values generated by `Nonpareil.curve()`_.

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

    source('utils/Nonpareil.R');
    samples <- read.table('samples.txt', sep='\t', h=T);
    attach(samples);
    np <- Nonpareil.curve.batch(File, 50, r=R, g=G, b=B, libnames=Name, modelOnly=TRUE);
    Nonpareil.legend('bottomright');
    detach(samples);


.. _R: http://www.r-project.org/

