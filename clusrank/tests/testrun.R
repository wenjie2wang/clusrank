## Please prepare a full list of test objects and test them one by one
## The current code does not even run.
library(clsrank)
data(crsd)
cluswilcox.test(z, cluster = id, data = crsd)
data(crsdUnb)
cluswilcox.test(z, cluster = id, data = crsdUnb)
data(crd)
cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
data(crdStr)
cluswilcox.test(z ~ group(group) + cluster(id) + stratum(stratum), data = crdStr)
