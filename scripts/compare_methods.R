# Functions for fitting other peoples' models
library(nnls)
library(provoc)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())


truth <- list(
    Deltacron = data.frame(
        variants = c("BA.1", "BA.2", "B.1.1.529", "B.1.617.2", "B.1.617.2+K417N"),
        probs = c(0.3, 0.2, 0.1, 0.35, 0.05)
    ),
    Delta_Plus_Plus = data.frame(
        variants = c("B.1.617.2", "B.1.617.2+K417N"),
        probs = c(0.97, 0.3)
    ),
    Too_Close = data.frame(
        variants = c("B.1.617.2", "B.1.617.2+K417N", "A.23.1", "A.23.1+E484K", 
            "AY.4.2", "AY.4"),
        probs = c(0.5, 0.03, 0.2, 0.02, 0.2, 0.05)
    ),
    Total_Below_One = data.frame(
        variants = c("BA.1", "BA.2", "B.1.1.529"),
        probs = c(0.2, 0.6, 0.2)
    ),
    One_Two_Punch = data.frame(
        variants = c("BA.1", "BA.2"),
        probs = c(0.2, 0.8)
    ),
    EZPZ = data.frame(
        variants = c("B.1.1.529", "B.1.617.2"),
        probs = c(0.5, 0.5)
    )
)

for (scenario in names(truth)) {
    print(scenario)
    varmat <- astronomize()
    varmat <- varmat[rownames(varmat) %in% truth[[scenario]]$variants,, drop = FALSE]
    varmat <- varmat[, apply(varmat, 2, sum) > 0]
    rownames(varmat) <- gsub("\\+", "_", rownames(varmat))
    varmat <- varmat[order(rownames(varmat)), ]
    true_vals <- data.frame(
        variant = row.names(varmat),
        prob = truth[[scenario]]$probs)
    rel_counts <- round(true_vals$prob * 1000)
    if(scenario == "Total_Below_One") {
        varmat2 <- varmat[2:3,, drop = FALSE]
        rel_counts2 <- rel_counts[2:3]
    } else {
        varmat2 <- varmat
        rel_counts2 <- rel_counts
    }

    for(i in 1:100) {
        coco <- simulate_coco(varmat, rel_counts = rel_counts, verbose = FALSE)
        fused <- fuse(coco, varmat2, verbose = FALSE)
        fused <- filter(fused, coverage > 0)
        # Fuse to ensure same mutation list with correct order
        cocovar <- provoc:::fission(fused)
        coco <- cocovar$coco
        varmat2 <- cocovar$varmat
        all_res_tmp <- bind_rows(
            alcov(coco, varmat2),
            optim_methods(coco, varmat2),
            avg_freq(fused, method = "Simple Avg"),
            avg_freq(fused, method = "Simple Med"),
            avg_freq(fused, method = "Binomial_GLM"),
            avg_freq(fused, method = "Binomial_GLM_Quasi")
        )
        all_res_tmp$iter <- i
        if(i == 1) {
            all_res <- all_res_tmp
        } else {
            all_res <- bind_rows(all_res, all_res_tmp)
        }
    }

    true_vals2 <- bind_rows(true_vals, 
        data.frame(variant = "total", 
            prob = sum(rel_counts2) / sum(rel_counts)))
    all_res1 <- all_res %>%
        group_by(method, iter) %>%
        summarise(variant = "total",
            rho = sum(rho, na.rm = TRUE), 
            .groups = "drop") %>%
        bind_rows(all_res)
    all_res2 <- left_join(all_res1, true_vals2, by = "variant")
    all_res2$scenario <- scenario
    if(scenario == names(truth)[1]) {
        all_res3 <- all_res2
    } else {
        all_res3 <- bind_rows(all_res3, all_res2)
    }
}

ggplot(all_res3, aes(x = method, y = rho, fill = method)) +
    geom_violin(draw_quantiles = c(0.045, 0.5, 0.954)) +
    facet_grid(scenario ~ variant, scales = "free") +
    geom_hline(aes(yintercept = prob, group = variant)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(all_res3, aes(x = method, y = rho - prob, fill = method)) +
    geom_violin(draw_quantiles = c(0.045, 0.5, 0.954)) +
    facet_grid(scenario ~ variant, scales = "free") +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

all_res3 %>%
    filter(!is.na(rho), !is.na(method), !is.na(scenario), !is.na(prob)) %>%
    group_by(method, scenario, variant) %>%
    summarise(abs_bias = abs(rho - prob), mse = mean((rho - prob)^2), var = var(rho), .groups = "drop") %>%
    group_by(method, scenario) %>%
    summarise(abs_bias = mean(abs_bias), mse = mean(mse), var = mean(var)) %>%
    #tidyr::pivot_longer(cols = c(abs_bias, var, mse)) %>%
    ggplot() +
        aes(x = method, y = mse, fill = method) +
        geom_col(position = "dodge", colour = "white") +
        facet_wrap( ~ scenario , scales = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

all_res3 %>%
    filter(variant != "total") %>%
    group_by(scenario, iter, method) %>%
    summarise(time = sum(time, na.rm = TRUE), .groups = "drop") %>%
    ggplot() +
        aes(x = method, y = time, fill = method) %>%
        geom_violin(scale = "width") +
        facet_wrap(~ scenario) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        labs(y = "Time (Seconds)")
