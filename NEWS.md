# clusrank 1.0-4

## Minor changes

* Added continuity correction to avoid zero p-value in permutation tests.
* Enhanced `clusWilcox.test()` to allow a `language` object.
  [#1](https://github.com/wenjie2wang/clusrank/issues/1)


# clusrank 1.0-3

## Minor changes

* Enhanced `clusWilcox.test()` to allow a `language` object.
  [#1](https://github.com/wenjie2wang/clusrank/issues/1)


# clusrank 1.0-2

## Bug fixes

* Fixed the tests for a nonzero location `mu` by removing a
  duplicated subtraction of `mu` from `x`.  Thank Iaroslav Pashentsev for the
  bug report.


# clusrank 1.0-1

## Bug fixes

* Fixed an internal function for encoding.  Thank Dipankar Bandyopadhyay for the
  bug report.
