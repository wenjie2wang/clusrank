\name{sedlab}
\alias{sedlab}
\docType{data}
\title{
Seddon Lab-CARMS-Wilcoxon data
}
\description{
Data from Seddon Lab. 
}
\usage{data(sedlab)}
\format{
  A data frame with 748 observations on the following 7 variables.
  \describe{
    \item{\code{ID}}{subject identifier}
    \item{\code{Eye}}{OD (right eye), OS(left eye)}
    \item{\code{Variant}}{\enumerate{
    \item no R1210C variant
    \item R1210C variant
    }
    }
    \item{\code{CARMS}}{Patient's last CARMS grade
    \enumerate{
     \item no AMD
     \item early AMD
     \item intermediate AMD
     \item geographic atroph(advanced dry)
     \item neovascular disease (advanced wet)  
    }
                                }
     \item{\code{Age_group}}{
     \enumerate{
     \item < 70 years
     \item 70 to 79.9 years 
     \item >= 80 years
     } 

     }
     \item{\code{Sex}}{
     \enumerate{
     \item male
     \item female
     }
     }
     \item{\code{Agesex}}{
     \enumerate{
     \item agegroup = 1; sex = 1
     \item agegroup = 2; sex = 1
     \item agegroup = 3; sex = 1
     \item agegroup = 1; sex = 2
     \item agegroup = 2; sex = 2
     \item agegroup = 3; sex = 2
     }

     }
                           
  }
}

\note{
CARMS grades were assessed separately for the two advanced stages (4 and 5):	
1. CARMS 1,2,3, and 4 was assessed, 	
2. CARMS 1,2,3, and 5 was assessed 	
}

\source{
Seddon Lab 
}

\examples{
data(sedlab)
cluswilcox(CARMS ~ Variant, sedlab, id = ID, stratum = Agesex)
}
\keyword{datasets}
\references{
Ferrara, Daniela and Seddon, Johanna M. (2015)
  Phenotypic Characterization of Complement Factor H R1210C Rare Genetic Variant in Age-Related Macular Degeneration, \emph{JAMA Ophthalmology} \bold{133},  785-791.
}