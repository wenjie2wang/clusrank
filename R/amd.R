################################################################################
##
## clusrank: Wilcoxon Rank Tests for Clustered Data
## Copyright (C) 2015-2022  Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
##
## This file is part of the R package clusrank.
##
## The R package clusrank is free software: You can redistribute it and/or
## modify it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or any later
## version (at your option). See the GNU General Public License at
## <https://www.gnu.org/licenses/> for details.
##
## The R package clusrank is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################
#' CARMS scores
#'
#'A data set from a research on complement factor H R1210C rare variant
#'and its associated phenotype.
#'This data set contains Clinical
#'Age-Related Maculopathy Staging (CARMS) scores
#'from a total of 143 patients (283 eyes), including
#'62 patients with the rare variant,
#'The data is from the lab of Dr. Johanna M. Seddon,
#'
#'
#'
#'
#' @format A data frame with 283 rows and 7 variables.
#'
#' \itemize{
#'   \item  ID patient identifier
#' \item  Eye OD (right eye), OS(left eye)
#' \item  Variant 1: No R1210C Variant; 2: R1210C Variant
#' \item CARMS Patient's last CARMS grade, related to age-related
#'     macular degeneration (AMD). 1: no AMD; 2: early AMD; 3:
#'     intermediate AMD;  4: geographic atroph(advanced dry); 5:
#'     neovascular disease (advanced wet)
#' \item Age_group 1: < 70 years;  2: 70 to 79.9 years; 3: >= 80 years
#' \item Sex 1: male; 2: female
#' \item AgeSex 1: agegroup = 1, sex = 1;  2: agegroup = 2, sex = 1;
#'     3: agegroup = 3, sex = 1; 4: agegroup = 1, sex = 2; 5: agegroup
#'     = 2, sex = 2; 6: agegroup = 3, sex = 2;
#' }
#'
#' @name amd
#' @docType data
#' @note CARMS grades were assessed separately for
#' the two advanced stages (4 and 5):
#' 1. CARMS 1,2,3, and 4 was assessed;
#' 2. CARMS 1,2,3, and 5 was assessed
#' @source The data came from Seddon's lab.
#' @references  Seddon JM, Sharma S, Adelman RA (2006)
#' \emph{Evaluation of the Clinical Age-related
#'  Maculopathy Staging System.}
#' Ophthalmology, \bold{113}, 260-266.
#' @references Ferrara D, Seddon JM (2015)
#' \emph{ Phenotypic characterization of complement
#' factor H R1210C rare genetic variant in
#' age-related macular degeneration}
#' JAMA Ophthalmol, 2015 Apr 16.
#' doi:10.1001/jamaophthalmol.2015.0814.
#' @examples
#' data(amd)
#' clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd,
#'                subset = CARMS %in% c(1, 2, 3, 4), method = "rgl", alternative = "two")
#' clusWilcox.test(CARMS ~ Variant + cluster(ID), data = amd,
#'                subset = CARMS %in% c(1, 2, 3, 4), method = "ds", alternative = "two")
#' clusWilcox.test(CARMS ~ Variant + cluster(ID) + stratum(AgeSex), data = amd,
#'                subset = CARMS %in% c(1, 2, 3, 4), alternative = "two")
#' @keywords datasets

NULL
