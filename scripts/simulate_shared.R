
library(nnls)
library(provoc)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(here)

source(here("scripts", "alt_methods.R"))

for (i in 0:19) {
    X1 <- c(rep(1, 20), rep(0, 20 - i))
    X2 <- c(rep(0, 20 - i), rep(1, 20))
    X3 <- rbinom(length(X1), prob = 0.5, size = 1)
    varmat <- do.call(rbind, list(X1, X2, X3))
    rownames(varmat) <- paste0("X", 1:nrow(varmat))
    colnames(varmat) <- paste0("M", 1:ncol(varmat))

    for(ii in 1:100) {
        print(ii)
        coco <- simulate_coco(varmat, rel_counts = c(200, 400, 50), verbose = FALSE)
        fused <- fuse(coco, varmat, verbose = FALSE)
        print("alcov")
            alcov(coco, varmat)
            print("optim")
            optim_methods(coco, varmat, method = c("freyja", "provoc"))
            print("avg")
            avg_freq(fused, method = "Simple Avg")
            avg_freq(fused, method = "Simple Med")
        
    }
}




