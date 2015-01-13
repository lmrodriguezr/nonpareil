Updates history
===============
What's new in 2.4
-----------------
* Subsampling now defaults to logarithmic: Previous versions subsampled linearly, since 2.4 the -d
  option defaults to 0.7. Options -m, -M, & -i still exist but they are ignored unless -d is 0.
* Nonpareil diversity: A logarithmic value of diversity is now reported, indicating the horizontal
  position of the Nonpareil curves. This value is the estimated mode of the fitted Gamma CDF, and
  cannot be calculated if model fitting fails.
* Experimental dataset comparisons: The -q option allows the paired comparison of datasets. The
  expectation is that the difference between the diversity of a sample by itself (without using -q)
  and the diversity of the same sample queried with a second sample (with -q) represents the distance
  between them. Note that this method experimental and is not symetric (hence not a real distance).
  Usually, the average of the distances in both ways can be used for clustering of samples.

What's new in 2.3
-----------------
* Updates history.
* Nonpareil R package: No changes in the functions were introduced, but they are now available
  as an R package complete with documentation. See :doc:`curves`.
* Make's ``install`` target: This includes the installation of binaries, the manual page, and the
  new Nonpareil R package.
* Simplified help: Only mandatory and commonly used arguments are now displayed with ``nonpareil -h``.
  Complete documentation is maintained both in the manual page (``man nonpareil``) and the
  `online documentation`_. In order to maintain a centralized documentation, the complete help messages
  for the Nonpareil R package are now self-contained, and were removed from the `online documentation`_.

.. _Online documentation: http://nonpareil.readthedocs.org/
