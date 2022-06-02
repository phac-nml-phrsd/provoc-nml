library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(provoc)
library(nnls)

params <- read.csv("data/clean/variant_params.csv")
N <- 50
total_seqs <- 1000
weeks <- unique(params$week)

for(week in weeks) {
    print(Sys.time())
    print(week)
    these_pars <- params[params$week == week,]
    pango_vars <- astronomize()
    delta_vars <- astronomize(path = "../pango-designation/lineage_constellations")
    pango_vars <- pango_vars[!rownames(pango_vars) %in% delta_vars,]
    varmat <- as.matrix(bind_rows(as.data.frame(pango_vars), as.data.frame(delta_vars)))
    varmat[is.na(varmat)] <- 0
    varmat_data <- varmat[row.names(varmat) %in% these_pars$lineage,]

    Other_Delta <- grep("AY", row.names(delta_vars), value = TRUE)
    Other_Delta <- Other_Delta[!Other_Delta %in% rownames(varmat_data)]
    varmat_true <- varmat[rownames(varmat) %in% rownames(varmat_data) | 
        rownames(varmat) %in% Other_Delta,]

    for(replicate in 1:N) {
        # Add in "Other Delta" completelty at random
        relative <- data.frame(lineage = rownames(varmat_true)) %>%
            left_join(these_pars[, c("lineage", "percent")], by = "lineage")
        to_fill <- sum(is.na(relative$percent))
        percs <- runif(sum(is.na(relative$percent)))
        percs <- these_pars$percent[these_pars$lineage == "Other Delta"] * percs / sum(percs)
        relative$percent[is.na(relative$percent)] <- percs

        coco <- simulate_coco(varmat_true, rel_counts = total_seqs * relative$percent, verbose = FALSE)
        fused <- fuse(coco, varmat_data, verbose = FALSE)
        cocovar <- provoc:::fission(fused)
        coco <- cocovar$coco
        varmat <- cocovar$varmat
        all_res_tmp <- bind_rows(
            alcov(coco, varmat, method = "AlCoV-LM"),
            alcov(coco, varmat, method = "AlCoV-Robust"),
            alcov(coco, varmat, method = "AlCoV-NNLS"),
            freyja(coco, varmat),
            avg_freq(fused, method = "Simple Avg"),
            #avg_freq(fused, method = "binomial"),
            #avg_freq(fused, method = "quasibinomial"),
            provoc_optim2(coco, varmat)
        )
        all_res_tmp$iter <- replicate
        all_res_tmp$week <- week
        all_res_tmp$week <- as.character(all_res_tmp$week)
        all_res_tmp2 <- all_res_tmp %>%
            group_by(method, iter, week) %>%
            summarise(rho = sum(rho), variant = rep("total", n()), .groups = "drop") %>%
            bind_rows(all_res_tmp)
        if(replicate == 1 & week == weeks[1]) {
            all_res <- all_res_tmp2
        } else {
            all_res <- bind_rows(all_res, all_res_tmp2)
        }
    }
}
params1 <- select(params, variant = lineage, true_perc = percent, week)
params2 <- params1 %>%
    filter(variant %in% rownames(varmat_data)) %>%
    group_by(week) %>%
    summarise(variant = "total", true_perc = sum(true_perc), .groups = "drop") %>%
    bind_rows(params1)
all_res2 <- left_join(all_res, params2, by = c("variant", "week"))

ggplot(all_res2) +
    aes(x = method, y = rho, colour = method) + 
    geom_hline(yintercept = c(0,1), colour = "darkgrey", linetype = "dashed") + 
    geom_violin() + geom_hline(mapping = aes(yintercept = true_perc)) + 
    facet_grid(week ~ variant) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
