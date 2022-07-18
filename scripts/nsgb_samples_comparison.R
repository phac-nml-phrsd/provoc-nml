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

truths <- list(
    deltacron100 = c("B.1.1.529" = 100, "B.1.617.2" = 100, "AY.25" = 100),
    deltacron1000 = c("B.1.1.529" = 1000, "B.1.617.2" = 1000, "AY.25" = 1000),
    deltacron10000 = c("B.1.1.529" = 10000, "B.1.617.2" = 10000, "AY.25" = 10000),
    oh_my_cron = c("B.1.1.529" = 500, "BA.1" = 1000, "BA.2" = 100),
    dealt_a_delta = c("B.1.617.2" = 100, "AY.4" = 1000, "AY.4.2" = 100, 
        "AY.25" = 500),
    new_delta = c("B.1.617.2" = 100, "BA.1" = 1000, "BA.2" = 1000),
    trace_amounts = c("B.1.617.2" = 50, "AY.4" = 50, "AY.4.2" = 50,
        "AY.25" = 50, "BA.1" = 2000, "BA.2" = 1000)
    )


# Setting up various variant matrices
constellations <- astronomize()
constellations <- constellations[rownames(constellations) %in% 
    c("AY.4", "AY.4.2", "B.1.1.529", "B.1.617.2", "BA.1", "BA.2"),]
constellations <- constellations[, colSums(constellations) > 0]
colnames(constellations)[which(colnames(constellations) == "+2205.GAGCCAGAA")] <- "ins:22205:9"
colnames(constellations)[which(colnames(constellations) == "+8262.AACA")] <- "ins:28262:4"
colnames(constellations)[which(colnames(constellations) == "28271-")] <- "del:28271:1"
# Spike protein starts at position 21562
# Position reported as index of amino acids, hence 3*246
colnames(constellations)[which(colnames(constellations) == "aa:S:RSYLTPG246-")] <- paste0("del:", 21562 + 3*246, ":21")
colnames(constellations)[which(colnames(constellations) == "aa:S:Y144-")] <- paste0("del:", 21562 + 3*144, ":1")
colnames(constellations)[which(colnames(constellations) == "aa:S:HV69-")] <- paste0("del:", 21562 + 3*69, ":2")
# ORF 1a starts at 265
colnames(constellations)[which(colnames(constellations) == "aa:orf1a:SGF3675-")] <- paste0("del:", 265 + 3*3675 - 1, ":9")

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
df2 <- df[, -1]
df2 <- df2[, -which(colnames(df2) %in% c("AY.25.1"))]

empirical <- apply(df2, 1, function(x) (x > 0.85) * (sum(x < 0.05) >= length(x) - 3))

rowSums(empirical)
rowSums(constellations)
rowSums(gromstollations)

varmat_list <- list(constellations, gromstollations, empirical)


# Sampling (currently very slow.)

sample_from_nsgb <- function(variants) {
    for(i in seq_along(variants)) {
        f <- read.csv(here("output", paste0(names(variants[i]), "_samplemat.csv")))
        ds <- which(names(f) == "date_submitted")
        count <- colSums(f[sample(1:nrow(f), size = variants[i], replace = TRUE), -c(1:ds)])
        df <- data.frame(count = count, mutation = names(count))

        if(i == 1) {
            res <- df
        } else {
            res <- full_join(res, df, by = "mutation")
            res[is.na(res)] <- 0
            res$count <- res$count.x + res$count.y
            res <- res[, c("mutation", "count")]
        }
    }

    res$coverage <- sum(variants)
    degredation <- runif(nrow(res), 0, 1)
    res$count <- round(res$count * degredation, 0)
    res$coverage <- round(res$coverage * degredation, 0)
    res <- res[res$coverage > 0, ]
    res$frequency <- round(res$count / res$coverage, 5)
    res
}





# Simulation and Fitting
t_total <- Sys.time()
for(rep in 1:500) {
    print(rep)
    for(sitch in seq_along(truths)) {
        print(sitch)
        situation <- truths[[sitch]]
        true_vals <- data.frame(variant = names(situation),
            rho = situation / sum(situation),
            se = NA, method = "Truth")

        for(varmat_type in 1:3) {
            t_sim0 <- Sys.time()
            coco <- sample_from_nsgb(situation)
            t_sim <- difftime(Sys.time(), t_sim0, units = "secs")
            print(t_sim)
            true_vals$time <- as.numeric(t_sim)

            coco$mutation <- gsub(".", ":", coco$mutation, fixed = TRUE)
            varmat <- varmat_list[[varmat_type]]
            fused <- fuse(coco, varmat, verbose = FALSE)
            fissed <- provoc:::fission(fused)
            coco <- fissed$coco
            varmat <- fissed$varmat


            all_res_tmp <- bind_rows(
                alcov(coco, varmat, method = c("AlCoV-LM", "AlCoV-Robust")),
                optim_methods(coco, varmat, method = c("freyja", "provoc", "squared")),
                avg_freq(fused, method = "Simple Avg"),
                avg_freq(fused, method = "Simple Med")
            )

            all_res_tmp <- bind_rows(all_res_tmp, true_vals)
            all_res_tmp$situation <- names(truths)[sitch]

            if(!exists("all_res")) {
                all_res <- all_res_tmp
            } else {
                all_res <- bind_rows(all_res, all_res_tmp)
            }

        }
    }
}
print(difftime(Sys.time(), t_total))

head(all_res)
all_res2 <- all_res
rm(all_res)

ggplot(all_res2) +
    aes(x = variant, y = rho, colour = method) +
    geom_violin() +
    facet_wrap(~ situation) +
    coord_flip()
