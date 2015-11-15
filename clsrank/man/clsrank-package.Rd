\name{clsrank-package}
\alias{clsrank-package}
\alias{clsrank}
\docType{package}
\title{
Wilcoxon rank sum test for clustered data
}
\description{
Incorporating clustering effects for the Wilcoxon Rank Sum 
Test for stratified balanced or unbalanced designs. Can be applied to
both stratified data and non-stratified data.
}
\details{
\tabular{ll}{
Package: \tab clsrank\cr
Type: \tab Package\cr
Version: \tab 0.0-4\cr
Date: \tab 2015-5-8\cr
License: \tab GPL (>=3)
\cr
}
Included is a function, cluswilcox, the arguments of this function is a formula (score ~ group), id, stratum(if any), data set name. Also two data sets are included: crd is a clustered non-stratified data, crd.str is a clustered stratified data.
}
\author{
Yujing Jiang <yujing.jiang@uconn.edu> [aut, cre], Jun Yan <yan.jun@uconn.edu> [ctb] 
Maintainer: Yujing.Jiang <yujing.jiang@uconn.edu>
}
\references{
Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2003) \emph{Incorporation of Clustering Effects for the Wilcoxon Rank Sum Test: A Large-Sample Approach}. Biometrics, \bold{59}, 1089-1098.
Bernard Rosner, Robert J. Glynn, Mei-Ling Ting Lee(2006) 
\emph{The Wilcoxon Signed Rank Test for Paired Comparisons of
 Clustered Data.} Biometrics, \bold{62}, 185-192.
}

\keyword{ package }
\examples{
data(crsd)
clusignrank(z, id, data = crsd)
data(crsd.unb)
clusignrank(z, id, data = crsd.unb)

data(crd)
cluswilcox(z ~ group, data = crd, id = id)
data(crd.str)
cluswilcox(z ~ group, data = crd.str, id = id, stratum = stratum)

}
