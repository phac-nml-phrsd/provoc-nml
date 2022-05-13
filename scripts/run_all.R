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

force_refit <- FALSE

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
varmat_type <- varmat_types[4]

for(varmat_type in varmat_types[1:4]) { # varmat_from_data causes errors - probably too big or with too many minor mutations

    handle <- paste0("results/", varmat_type, ".RDS")
    if(file.exists(handle) & !force_refit) {
        res <- readRDS(handle)
        res1 <- res$optim
        res2 <- res$runjags
    } else {
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
        } else if (varmat_type == "varmat_from_data") {
            varmat <- varmat_from_data(type = coco$type, pos = coco$pos, alt = coco$alt, max_n = 30, mutation_format = "aa")
        } else {
            stop("No valid varmat type specified.")
        }

        fused <- fuse(coco, varmat)
        dim(coco); dim(fused)

        res1 <- provoc(fused = fused, method = "optim")
        res1$method <- "optim"
        res2 <- provoc(fused = fused, method = "runjags", quiet = 0)
        res2$method <- "runjags"
        saveRDS(list(optim = res1, runjags = res2), 
            file = handle)
    }
    df <- bind_rows(res1, res2)

    pdf(file = paste0("results/", varmat_type, ".pdf"), width = 11, height = 8.5)

    # Optim Results
    gg_optim <- res1 %>% group_by(variant) %>%
        mutate(include = rep(!all(rho < 0.05), n())) %>%
        ungroup() %>% 
        filter(include) %>%
        ggplot(aes(x = ymd(date), y = rho, colour = location)) + 
            geom_point() + 
            facet_wrap(~ variant) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
            scale_x_date(breaks = sort(unique(ymd(res1$date)))) +
            labs(title = "Results from Optim")
    print(gg_optim)

    # RunJags Results
    gg_runjags <- res2 %>% group_by(variant) %>%
        mutate(include = rep(!all(rho < 0.05), n())) %>%
        ungroup() %>% 
        filter(include) %>%
        ggplot(aes(x = ymd(date), y = rho, ymin = ci_low, ymax = ci_high, colour = location)) + 
            geom_point() + geom_errorbar() +
            facet_wrap(~ variant) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
            scale_x_date(breaks = sort(unique(ymd(res1$date)))) +
            labs(title = "Results from RunJAGS")
    print(gg_runjags)

    # Both methods at once
    gg_both <- ggplot(df) + 
        aes(x = ymd(date), y = rho, shape = method, colour = location) +
        geom_point() + 
        geom_line(aes(group = paste(method, location))) +
        facet_wrap(~ variant) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        scale_x_date(breaks = sort(unique(ymd(res1$date)))) +
        labs(title = "Results Everywhere")
    print(gg_both)

    # Both methods - Stacked Bar Chart
    gg_stacked <- ggplot(res1) +
        aes(x = ymd(date), y = rho, fill = variant) +
        geom_col(colour = 1) +
        facet_wrap(~ method + location) +
        coord_flip()  +
        scale_x_date(breaks = sort(unique(ymd(res1$date)))) +
        labs(title = "Results Everywhere")
    print(gg_stacked)

    # Concordance between methods
    gg_conc <- full_join(res1, res2, by = c("date", "location", "sample", "variant")) %>%
        ggplot() +
            aes(x = rho.x, y = rho.y) +
            geom_point() + 
            geom_abline(intercept = 0, slope = 1) +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
            labs(x = "Optim", y = "RunJags") +
            facet_wrap(~ variant) +
            labs(title = "Concordance between methods")
    print(gg_conc)
    dev.off()

    convergence(res1)
    convergence(res2)
}
