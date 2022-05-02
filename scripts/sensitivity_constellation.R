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

library(parallel)

animal <- read.csv(here("data/clean", "nml.csv")) 
coco <- animal %>%
    mutate(location = ifelse(grepl("MMN", sample), 
        yes = "MMN", 
        no = ifelse(grepl("MMS", sample),
            yes = "MMS",
            no = "VLI")),
        mutation = label,
        coverage = alt_dp + ref_dp, # TODO: Is this correct?
        count = round(frequency * coverage, 0))

coco <- filter(coco, location == "VLI")

varmat <- astronomize(path = here("..", "/constellations"))


voi <- c("B.1.617.2", "B.1.617.2+K417N", "AY.4")

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

par_fun <- function(i) {
    iteration <- inu[i, 'iteration']
    nuisances <- inu[i, 'nuisance']
    oldvars <- rownames(varmat)[!rownames(varmat) %in% voi]
    newvars <- sample(oldvars, nuisances, FALSE)
    varmat2 <- varmat[rownames(varmat) %in% c(voi, newvars), ]

    fused2 <- fuse(coco, varmat2, verbose = FALSE)

    summarise(group_by(mutate(provoc(fused = fused2, method = "optim"), 
            variant = ifelse(variant %in% voi, variant, "Other"),
            iteration = rep(iteration, n()),
            nuisance_lineages = rep(nuisances, n())
        ),
        variant, sample, date, iteration, nuisance_lineages),
        rho = sum(rho))
}

inu <- expand.grid(iteration = 1:100, nuisance = c(5, 10, 15, 20))

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, c("varmat", "coco", "voi", "provoc", "fuse", 
    "summarise", "group_by", "mutate", "nuisances", "n", "inu"))
allsense <- bind_rows(parLapply(cl, 1:nrow(inu), par_fun))
stopCluster(cl)

ggplot(allsense) + 
    aes(x = variant, y = rho, fill = factor(nuisance_lineages)) +
    geom_violin(draw_quantiles = c(0.5)) + 
    facet_wrap(~ variant + date, scales = "free_x") +
    labs(fill = "Nuisances")
# The increase in other as we add lineages indicates more than just Delta!
