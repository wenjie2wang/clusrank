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
