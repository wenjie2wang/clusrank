

library(clusrank)
data(crsd)
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl")
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "ds")
data(crsdUnb)
clusWilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "rgl")
clusWilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "ds")
data(crd)
clusWilcox.test(z ~ group + cluster(id), data = crd)
data(crdStr)
clusWilcox.test(z ~ group + cluster(id) + stratum(stratum), data = crdStr)
