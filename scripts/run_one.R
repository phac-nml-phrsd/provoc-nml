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

varmat_types <- c("constellations", "varmat_from_variants-all_voc", "varmat_from_variants-omicron_delta", "varmat_from_data")
# Choose a variant matrix
varmat_type <- varmat_types[4]

# Choose a method
method <- c("optim", "runjags")[1]

# (Optional) Choose one sample to run.
coco <- filter(coco, sample == unique(coco$sample)[1])


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
    varmat <- varmat_from_data(type = coco$type, pos = coco$pos, alt = coco$alt, max_n = 30, mutation_format = "aa")
}

fused <- fuse(coco, varmat)
res <- provoc(fused = fused, method = method)

res %>% 
    # Uncomment to only show high probability variants
    #group_by(variant) %>% mutate(include = rep(!all(rho < 0.05), n())) %>% ungroup() %>%  filter(include) %>%
    ggplot(aes(x = ymd(date), y = rho, colour = location)) + 
        geom_point() + 
        facet_wrap(~ variant) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        scale_x_date(breaks = sort(unique(ymd(res$date)))) +
        labs(title = "Results from Optim")
