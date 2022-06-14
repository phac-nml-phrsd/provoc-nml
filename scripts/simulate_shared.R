
library(nnls)
library(provoc)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(here)

source(here("scripts", "alt_methods.R"))

for (i in 0:19) {
    print(i)
    X1 <- c(rep(1, 20), rep(0, 20 - i), rep(0, 5))
    X2 <- c(rep(0, 20 - i), rep(1, 20), rep(0, 5))

    for(ii in 1:1000) {
        X3 <- rbinom(length(X1), prob = 0.5, size = 1)
        varmat <- do.call(rbind, list(X1, X2, X3))
        rownames(varmat) <- paste0("Variant", 1:nrow(varmat))
        colnames(varmat) <- paste0("M", 1:ncol(varmat))
        coco <- simulate_coco(varmat, rel_counts = c(420, 120, 60), verbose = FALSE)
        fused <- fuse(coco, varmat, verbose = FALSE)
        res <- bind_rows(
            alcov(coco, varmat, method = c("AlCoV-LM", "AlCoV-Robust", "AlCoV-Binom")),
            optim_methods(coco, varmat),
            avg_freq(fused, method = "Simple Avg"),
            avg_freq(fused, method = "Simple Med")
            )
        res$overlap <- i
        if(ii == 1) {
            res_tmp <- res
        } else {
            res_tmp <- bind_rows(res_tmp, res)
        }
    }

    if(i == 0) {
        all_res <- res_tmp
    } else {
        all_res <- bind_rows(all_res, res_tmp)
    }
}

head(all_res)

ggplot(all_res) +
    aes(x = factor(overlap, ordered = TRUE),
        y = rho, fill = method) + 
    geom_violin(position = position_dodge()) +
    facet_grid(method ~ variant)

all_res_stats <- all_res %>%
    group_by(variant, method, overlap) %>%
    summarise(m = mean(rho, na.rm = TRUE), 
        s = sd(rho, na.rm = TRUE),
        t = mean(time))

ggplot(all_res_stats) +
    aes(x = overlap, y = t, colour = method) +
    geom_point() + geom_line() +
    facet_wrap(~ variant)
ggplot(filter(all_res_stats, method %in% c("freyja", "AlCoV-LM", "Simple Avg", "Simple Med", "provoc"))) +
    aes(x = overlap, y = m, colour = method) +
    geom_point() + geom_line() +
    facet_wrap(~ variant)

ggplot(filter(all_res_stats, method %in% c("freyja", "AlCoV-LM", "Simple Avg", "Simple Med", "provoc"))) +
    aes(x = overlap, y = s, colour = method) +
    geom_point() + geom_line() +
    facet_wrap(~ variant) +
    labs(x = "Number Overlapping Mutations",
        y = "Standard Deviation of Estimates")
