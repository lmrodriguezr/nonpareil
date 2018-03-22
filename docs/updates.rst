Updates history
===============

What's new in 3.3
-----------------
* Release used in the manuscript "Nonpareil 3: Fast Estimation of Metagenomic
  Coverage and Sequence Diversity".

What's new in 3.2
-----------------
* Minor bug fixes and documentation updates.

What's new in 3.1
-----------------
* Implemented additional controls for the robustness of the Nonpareil sequence
  diversity index (Nd).
* Revamped R code, now with Object-oriented structure and uniform R-style
  options. The new package also includes test data and examples and is now
  available at `CRAN <https://CRAN.R-project.org/package=Nonpareil>`.

What's new in 3.0
-----------------
* New k-mer kernel: A much faster kernel is now provided as an alternative to
  the alignment-based kernel in previous versions.

What's new in 2.4
-----------------
* Subsampling now defaults to logarithmic: Previous versions subsampled
  linearly, since 2.4 the -d option defaults to 0.7. Options -m, -M, & -i still
  exist but they are ignored unless -d is 0.
* Nonpareil diversity: A logarithmic value of diversity is now reported,
  indicating the horizontal position of the Nonpareil curves. This value is the
  estimated mode of the fitted Gamma CDF, and cannot be calculated if model
  fitting fails.
* Experimental dataset comparisons: The -q option allows the paired comparison
  of datasets. The expectation is that the difference between the diversity of a
  sample by itself (without using -q) and the diversity of the same sample
  queried with a second sample (with -q) represents the distance between them.
  Note that this method experimental and is not symetric (hence not a real
  distance). Usually, the average of the distances in both ways can be used for
  clustering of samples.

What's new in 2.3
-----------------
* Updates history.
* Nonpareil R package: No changes in the functions were introduced, but they are
  now available as an R package complete with documentation. See :doc:`curves`.
* Make's ``install`` target: This includes the installation of binaries, the
  manual page, and the new Nonpareil R package.
* Simplified help: Only mandatory and commonly used arguments are now displayed
  with ``nonpareil -h``. Complete documentation is maintained both in the manual
  page (``man nonpareil``) and the `online documentation`_. In order to maintain
  a centralized documentation, the complete help messages for the Nonpareil R
  package are now self-contained, and were removed from the
  `online documentation`_.

.. _Online documentation: http://nonpareil.readthedocs.org/
