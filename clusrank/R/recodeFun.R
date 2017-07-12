################################################################################
##
##   R package clusrank by Mei-Ling Ting Lee, Jun Yan, and Yujing Jiang
##   Copyright (C) 2015-2017
##
##   This file is part of the R package clusrank.
##
##   The R package clusrank is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package clusrank is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package clusrank. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
recoderFunc <- function(data, oldvalue, newvalue) {
    
    ## convert any factors to characters    
    if (is.factor(data))     data     <- as.character(data)
    if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
    if (is.factor(newvalue)) newvalue <- as.character(newvalue)
    
    ## create the return vector
    
    newvec <- data
    
    ## put recoded values into the correct position in the return vector
    
    for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
    
    newvec
    
}

