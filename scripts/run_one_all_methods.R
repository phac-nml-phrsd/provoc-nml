suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    theme_set(theme_bw())
    library(provoc)
    library(lubridate)
    library(here)
    library(rjags)
    library(runjags)
    library(nnls)
})

source(here("scripts", "alt_methods.R"))

varmat_type <- c("constellations", "varmat_from_variants-all_voc", 
    "varmat_from_variants-omicron_delta", "varmat_from_data")[1]

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
coco <- filter(coco_tmp, sample == unique(coco_tmp$sample)[4])
cover <- read.csv(
    paste0("data/coverage/", unique(coco_tmp$sample)[4], 
        ".ivar_trim.sorted.bam_depth.tsv"), 
    header = FALSE, sep = "\t")
names(cover) <- c("specimen", "position", "coverage")

varmat <- astronomize(path = "../constellations/constellations/definitions")
rownames(varmat) <- gsub("(\\+)|-", "_", rownames(varmat))

coco2 <- add_coverage(coco, cover, colnames(varmat))

fused <- fuse(coco2, varmat)

cocovar <- provoc:::fission(fused)
coco <- cocovar$coco
varmat2 <- cocovar$varmat
all_res_tmp <- bind_rows(
    alcov(coco, varmat2, method = c("AlCoV-LM", "AlCoV-Robust", "AlCoV-Binom")),
    optim_methods(coco, varmat2)
    #avg_freq(fused, method = "Simple Avg"),
    #avg_freq(fused, method = "Simple Med"),
    #avg_freq(fused, method = "Binomial_GLM"),
    #avg_freq(fused, method = "Binomial_GLM_Quasi")
)

# DANGER! AlCoV has values below 0, but I'm not sure if this is how they deal with them,
all_res_tmp$rho[all_res_tmp$rho < 0] <- 0

ggplot(all_res_tmp) +
    aes(x = method, y = rho, colour = method) +
    geom_hline(yintercept = 0) +
    geom_point() +
    scale_colour_brewer(palette = 2, type = "qual") +
    facet_wrap(~ variant) +
    coord_flip()

preds <- lapply(unique(all_res_tmp$method), function(x) {
    rho <- all_res_tmp$rho[all_res_tmp$method == x]
    if(length(rho) == nrow(varmat2)) {
        pred <- rho %*% varmat2
        data.frame(mutation = colnames(varmat2), 
            pred_freq = as.numeric(pred),
            method = x)
    }
})

all_preds <- bind_rows(preds) %>%
    left_join(select(coco, mutation, count, coverage), by = "mutation") %>%
    mutate(pred_count = pred_freq * coverage)
ggplot(all_preds) +
    aes(x = pred_count, y = count, colour = method, shape = method) + 
    scale_shape_manual(values = 1:10) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~ method)

all_preds %>%
    group_by(method) %>%
    summarise(rmse = sqrt(mean((count - pred_count)^2, na.rm = TRUE)),
        rmse_freq = sqrt(mean((count/coverage - pred_freq)^2, na.rm = TRUE)),
        .groups = "drop") %>%
    tidyr::pivot_longer(-method) %>%
    ggplot() +
        aes(x = method, y = value, colour = method) +
        geom_point(size = 5) + 
        geom_segment(aes(yend = 0, xend = method)) +
        facet_wrap(~ name, scales = "free") +
        theme(legend.position = "none") +
        coord_flip()
