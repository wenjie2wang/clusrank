

library(clusrank)
data(crsd)
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl")
clusWilcox.test(z, cluster = id, data = crsd, paired = TRUE, method = "rgl", alternative = "greater", mu = 1)
clusWilcox.test(z ~ cluster(id), data = crsd, paired = TRUE, method = "ds")
clusWilcox.test(z ~ cluster(id), data = crsd, paired = TRUE, method = "rgl", exact = TRUE)
data(crsdUnb)
clusWilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "rgl")
clusWilcox.test(z, cluster = id, data = crsdUnb, paired = TRUE, method = "ds")
data(crd)
clusWilcox.test(z ~ group + cluster(id), data = crd)
crd1 <- crd[c(1:20, 141:160), ]
clusWilcox.test(z ~ group + cluster(id), data = crd1, method = "rgl", exact = TRUE)
data(crdStr)
clusWilcox.test(z ~ group + cluster(id) + stratum(stratum), data = crdStr)
clusWilcox.test(z ~ group + cluster(id) + stratum(stratum), data = crdStr, method = "ds")
