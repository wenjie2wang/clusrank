# clusrank

[![CRAN_Status_Badge][r-pkg-badge]][cran-url]
[![Build Status][gha-icon]][gha-url]


The R package clusrank provides functions for Wilcoxon rank sum test and
Wilcoxon signed rank test for clustered data.
See Jiang et. al (2020) for details.

## Installation

You may install the released version from [CRAN][cran-url].

```R
install.packages("clusrank")
```

## Development

The latest version of package is under development at [GitHub][github-url].  If
it is able to pass the automated package checks, one may install it by

```R
if (! require(remotes)) install.packages("remotes")
remotes::install_github("wenjie2wang/clusrank", upgrade = "never")
```

## References

- Jiang, Y., He, X., Lee, M. T., Rosner, B., & Yan, J. (2020). Wilcoxon
  rank-based tests for clustered data with R package clusrank. Journal of
  Statistical Software, 96(6), 1–26. http://dx.doi.org/10.18637/jss.v096.i06


## License

[GNU General Public License][gpl] (≥ 3)


[r-pkg-badge]: https://www.r-pkg.org/badges/version/clusrank
[cran-url]: https://CRAN.R-project.org/package=clusrank
[github-url]: https://github.com/wenjie2wang/clusrank
[gha-icon]: https://github.com/wenjie2wang/clusrank/workflows/R-CMD-check/badge.svg
[gha-url]: https://github.com/wenjie2wang/clusrank/actions
[gpl]: https://www.gnu.org/licenses/
