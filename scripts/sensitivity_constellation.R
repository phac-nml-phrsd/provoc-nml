# Sensitivity Analysis for inclusion of variants

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    theme_set(theme_bw())
    library(provoc)
    library(lubridate)
    library(here)
    library(rjags)
    library(runjags)
})

animal <- read.csv(here("data/clean", "nml.csv")) 
coco <- animal %>%
    mutate(location = ifelse(grepl("MMN", sample), 
        yes = "MMN", 
        no = ifelse(grepl("MMS", sample),
            yes = "MMS",
            no = "VLI")),
        mutation = label,
        coverage = alt_dp + ref_dp, # TODO: Is this correct for coverage?
        count = round(frequency * coverage, 0))

coco <- filter(coco, location == "VLI")

varmat <- astronomize(path = here("..", "/constellations"))


voi <- c("B.1.617.2", "B.1.617.2+K417N", "AY.4", "AY.4.2")

sens10 <- bind_rows(lapply(1:100, function(i) {
    oldvars <- rownames(varmat)[!rownames(varmat) %in% voi]
    newvars <- sample(oldvars, 10, FALSE)
    varmat2 <- varmat[rownames(varmat) %in% c(voi, newvars), ]
    fused2 <- fuse(coco, varmat2)
    res <- provoc(fused = fused2, method = "optim")
    res %>% 
        mutate(
            variant = ifelse(variant %in% voi, variant, "Other10"),
            iteration = i,
            nuisance_lineages = 10
        ) %>%
        group_by(variant, sample, date) %>%
        summarise(rho = sum(rho))
}))


ggplot(sens10) + 
    aes(x = date, y = rho) +
    geom_violin() + 
    facet_wrap(~ variant) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# TODO: scatterplots of lineages to show correspondence

t0 <- Sys.time()
for(nuisances in c(5, 10, 15, 20, 25)) {
    for(i in 1:100) {
        print(nuisances)
        print(i)
        t1 <- Sys.time()
        oldvars <- rownames(varmat)[!rownames(varmat) %in% voi]
        newvars <- sample(oldvars, nuisances, FALSE)
        varmat2 <- varmat[rownames(varmat) %in% c(voi, newvars), ]

        fused2 <- fuse(coco, varmat2, verbose = FALSE)

        res <- provoc(fused = fused2, method = "optim") %>% 
            mutate(
                variant = ifelse(variant %in% voi, variant, "Other"),
                iteration = rep(i, n()),
                nuisance_lineages = rep(nuisances, n())
            ) %>%
            group_by(variant, sample, date, iteration, nuisance_lineages) %>%
            summarise(rho = sum(rho))

        if(!exists("allsense")) {
            allsense <- res
        } else {
            allsense <- bind_rows(allsense, res)
        }
        print(Sys.time() - t1)
        print()
    }
}
Sys.time() - t0

ggplot(allsense) + 
    aes(x = variant, y = rho, fill = factor(nuisance_lineages)) +
    geom_violin() + 
    facet_wrap(~ date + variant, scales = "free_x")
