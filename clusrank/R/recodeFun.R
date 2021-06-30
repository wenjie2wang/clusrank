################################################################################
##
## clusrank: Wilcoxon Rank Tests for Clustered Data
## Copyright (C) 2015-2021  Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
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
recoderFunc <- function(data, oldvalue, newvalue) {

    ## convert any factors to characters
    if (is.factor(data))     data     <- as.character(data)
    if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
    if (is.factor(newvalue)) newvalue <- as.character(newvalue)

    ## create the return vector

    newvec <- rep(NA, length(data))

    ## put recoded values into the correct position in the return vector

    for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]

    newvec

}
