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

varmat_type <- c("constellations", "varmat_from_variants-all_voc", 
    "varmat_from_variants-omicron_delta", "varmat_from_data")[1]
method <- c("optim", "runjags")[1]

animal <- read.csv(here("data/clean", "nml.csv")) 
coco_tmp <- animal %>%
    mutate(location = ifelse(grepl("MMN", sample), 
        yes = "MMN", 
        no = ifelse(grepl("MMS", sample),
            yes = "MMS",
            no = "VLI")),
        mutation = label,
        coverage = total_dp, # TODO: Is this correct for coverage?
        count = round(frequency * coverage, 0))


# Choose one sample to run.
coco <- filter(coco_tmp, sample == unique(coco_tmp$sample)[3])
cover <- read.csv(
    paste0("data/coverage/", unique(coco_tmp$sample)[3], 
        ".ivar_trim.sorted.bam_depth.tsv"), 
    header = FALSE, sep = "\t")
names(cover) <- c("specimen", "position", "coverage")

cococo <- left_join(coco, cover, by = c("pos" = "position"))
plot(cococo$coverage.x ~ cococo$coverage.y,
    xlab = "Coverage from coverage file", 
    ylab = "TOTAL_DP")
abline(0, 1) # There's a problem here - ignoring for now.


if(varmat_type == "constellations") {
    varmat <- astronomize(path = here("..", "/constellations"))
} else if (varmat_type == "varmat_from_variants-all_voc") {
    all_voc <-  c("B.1.1.529", "BA.1", "BA.1.1", "BA.2", 
        "B.1.1.7", "P.1", "P.2", "P.3", 
        "B.1.617.2", "AY.1", "AY.2", "AY.4", "AY.4.2", 
        "C.37", "B.1.621", "B.1.351", "B.1.525", "B.1.526", "B.1.617.1")
    all_voc <- all_voc[all_voc %in% unique(mutations_by_lineage$lineage)]
    varmat <- varmat_from_variants(variants = all_voc, 
        mutation_format = "aa")
} else if (varmat_type == "varmat_from_variants-omicron_delta") {
    varmat <- varmat_from_variants(variants = c("B.1.1.529", "BA.1", "BA.1.1", "BA.2", 
        "B.1.617.2", "AY.1", "AY.2", "AY.4", "AY.4.2"), 
        mutation_format = "aa")
} else {
    varmat <- varmat_from_data(type = coco$type, pos = coco$pos, alt = coco$alt, max_n = 80, mutation_format = "aa", matches = 3)
}

coco2 <- add_coverage(coco, cover, colnames(varmat))

fused <- fuse(coco, varmat)
res <- provoc(fused = fused, method = method)

res %>% 
    ggplot(aes(x = variant, y = rho)) + 
        geom_point() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(title = "Results from Optim") +
        coord_flip()

var2 <- varmat[res$variant,]

pred_muts <- res$rho %*% var2

pred_compare <- inner_join(
    x = coco[, c("mutation", "count", "coverage", "frequency")],
    y = data.frame(mutation = colnames(pred_muts), 
        rho = as.numeric(pred_muts)),
    by = "mutation")

ggplot(pred_compare) +
    aes(x = frequency, y = rho) + 
    geom_point()
