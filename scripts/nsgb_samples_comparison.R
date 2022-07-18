# NSGB re-sampling and analysing with the methods.

suppressPackageStartupMessages({
    library(nnls)
    library(provoc)
    library(dplyr)
    library(ggplot2)
    theme_set(theme_bw())
    library(here)
})

source(here("scripts", "alt_methods.R"))

truth <- list(
    deltacron100 = c("B.1.1.529" = 100, "B.1.617.2" = 100, "AY.25" = 100),
    deltacron1000 = c("B.1.1.529" = 1000, "B.1.617.2" = 1000, "AY.25" = 1000),
    deltacron10000 = c("B.1.1.529" = 10000, "B.1.617.2" = 10000, "AY.25" = 10000),
    oh_my_cron = c("B.1.1.529" = 500, "BA.1" = 1000, "BA.2" = 100),
    dealt_a_delta = c("B.1.617.2" = 100, "AY.4" = 1000, "AY.4.2" = 100, "AY.24" = 500, 
        "AY.25" = 100, "AY.25.1" = 1000, "AY.43" = 500),
    new_delta = c("B.1.617.2" = 100, "BA.1" = 1000, "BA.2" = 1000),
    trace_amounts = c("B.1.617.2" = 50, "AY.4" = 50, "AY.4.2" = 50, "AY.24" = 50, 
        "AY.25" = 50, "AY.25.1" = 50, "AY.43" = 50, "BA.1" = 2000, "BA.2" = 1000)
    )


# Setting up various variant matrices
constellations <- astronomize()
gromstollations <- astronomize(path = here("data", "gromstollations"))

# empirical mutation lists based on what I've sampled
nsgb <- list.files(here("output"), pattern = "_samplemat.csv", 
    full.names = TRUE)
for(i in seq_along(nsgb)) {
    f <- read.csv(nsgb[i])
    ds <- which(colnames(f) == "date_submitted")
    empiric_perc <- apply(f[, -c(1:ds)], 2, mean)
    empiric <- data.frame(mutation = names(empiric_perc), perc = empiric_perc)

    # Parse filename for lineage name
    names(empiric)[2] <- gsub("_samplemat.csv", "", rev(strsplit(nsgb[i], "/")[[1]])[1]) 

    if(i == 1) {
        df <- empiric
    } else {
        df <- full_join(df, empiric, by = "mutation")
    }
}
df[is.na(df)] <- 0
rownames(df) <- df$mutation
df <- df[, -1]
df <- df[, -which(colnames(df) %in% c("AY.25.1", "AY.4.2"))]

empirical <- apply(df, 1, function(x) (x > 0.85) * (sum(x < 0.05) >= length(x) - 3))



