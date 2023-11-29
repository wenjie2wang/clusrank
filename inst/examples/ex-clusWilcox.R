library(clusrank)

## Clustered signed rank test using RGL method.
data(crsd)
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl")
## or
clusWilcox.test(z ~ cluster(id), data = crsd, paired = TRUE, method = "rgl")

## Clustered rank sum test using RGL method.
data(crd)
clusWilcox.test(z ~ group + cluster(id), data = crd, method = "rgl")
## or
clusWilcox.test(z, cluster = id, group = group, data = crd, method = "rgl")
