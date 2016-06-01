

library(clusrank)
data(crsd)
cluswilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl")
cluswilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "ds")
data(crsdUnb)
cluswilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "rgl")
cluswilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "ds")
data(crd)
cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
data(crdStr)
cluswilcox.test(z ~ group(group) + cluster(id) + stratum(stratum), data = crdStr)
