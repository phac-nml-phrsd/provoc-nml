library(jsonlite)

all_voc <-  c("B.1.1.529", "BA.1", "BA.1.1", "BA.2", 
        "B.1.1.7", "P.1", "P.2", "P.3", 
        "B.1.617.2", "AY.1", "AY.2", "AY.4", "AY.4.2", 
        "C.37", "B.1.621", "B.1.351", "B.1.525", "B.1.526", "B.1.617.1")

deltemp_aa <- read_json("https://lapis.cov-spectrum.org/open/v1/sample/aa-mutations?pangoLineage=B.1.617.2")
deltemp_nuc <- read_json("https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage=B.1.617.2")

foo_aa <- data.frame(t(sapply(deltemp_aa[[3]], as.character)))
foo_nuc <- data.frame(t(sapply(deltemp_nuc[[3]], as.character)))

dim(foo_aa)
dim(foo_nuc)

foo_aa[,1][order(foo_aa[,2])]
foo_nuc[,1][order(foo_nuc[,2])]
