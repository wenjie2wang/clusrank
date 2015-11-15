\name{crd}
\alias{crd}
\docType{data}
\title{
Clustered non-stratified data
}
\description{
Test non-stratified data include 3 variables: id, z, group. Have some missing values.
}
\usage{data("crd")}
\format{
  A data frame with 340 observations on the following 3 variables.
  \describe{
    \item{\code{id}}{id of cluster}
    \item{\code{z}}{score of observation}
    \item{\code{group}}{The group each observation comes from}
  }
}
\source{
URL: https://sites.google.com/a/channing.harvard.edu/bernardrosner/channing/incorporating-cluster-effects-for-the-wilcoxon-rank-sum-test-1/cluswilcox/sample_data
}

\examples{
data(crd)
}
\keyword{datasets}
