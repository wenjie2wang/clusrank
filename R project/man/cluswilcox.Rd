\name{cluswilcox}
\alias{cluswilcox}
\title{
Wilcoxon rank sum test for clustered data
}
\description{
Incorporating clustering effects for the Wilcoxon rank sum test for
stratified balanced or unbalanced designs. Can be applied to both
stratified data and non-stratified data.
}
\usage{
cluswilcox(formula, data = parent.frame(), subset = NULL,
           id = NULL, stratum = NULL, contrasts = NULL,
           na.action = na.omit, ...)
}
\arguments{
  \item{formula}{
    an object of class \code{"formula"} in the form of score ~ group,
    where score is the response variable of interest and group is the
    grouping factor.    
  }
  \item{data}{
    an optional data frame, list or environment (or object
    coercible by \code{as.data.frame} to a data frame) containing the
    variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)}, typically
    the environment from which this function is called.
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be used.
  }
  \item{id}{
    a numeric vector indicating the cluster ids of the scores. If not
    specified, each score has its own id, which reduces to non-clustered
    data.
  }
  \item{stratum}{
    a numeric vector indicating the stratum ids of the scores. If not
    specified, the scores are from the same cluster, which reduces to
    non-stratified data.
  }
  \item{contrasts}{
    the same as that in \code{lm}. Currently unused.
  }
  \item{na.action}{
    a function which indicates what should happen when the data contains
    NAs. The  default action is to omit them.
  }
  \item{...}{
    additional arguments, currently ignored.
  }
}
\value{
  a list with the following components
  \item{WC}{Clustered Wilcoxon ranksum statistic}
  \item{ExpWc}{Expected value clustered Wilcoxon ranksum sstatistic}
  \item{varWc}{Variance of clustered Wilcoxon ranksum statistic}
  \item{zc}{z statistic for clustered Wilcoxon ranksum statistic}
  \item{p.value}{P-value for clustered Wilcoxon ranksum z statistic}
}
\references{
  Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2003) \emph{Incorporation of Clustering Effects for the Wilcoxon Rank Sum Test: A Large-Sample Approach.} Biometrics, \bold{59}, 1089-1098.
}
\author{
  Yujing Jiang <yujing.jiang@uconn.edu>
}
\examples{
data(crd)
cluswilcox(z ~ group, data = crd, id = id)
data(crd.str)
cluswilcox(z ~ group, data = crd.str, id = id, stratum = stratum)
}
