##
## clusrank: Wilcoxon Rank Tests for Clustered Data
##
## Copyright (C) 2015-2024 Yujing Jiang, Mei-Ling Ting Lee, and Jun Yan
## Copyright (C) 2022-2024 Wenjie Wang
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


## returns the p-value for permutation tests
perm_pvalue <- function(w0, ws, alternative)
{
    ## continuity correction to prevent p-value from being exactly zero
    prob_le_w <- (sum(ws <= w0) + 0.5) / (length(ws) + 1)
    switch(
        alternative,
        two.sided = 2 * min(prob_le_w, 1 - prob_le_w),
        greater = 1 - prob_le_w,
        less = prob_le_w,
        stop("Unknown alternative.")
    )
}
