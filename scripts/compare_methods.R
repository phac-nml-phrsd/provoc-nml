# Functions for fitting other peoples' models
library(nnls)
library(provoc)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

alcov <- function(coco, varmat, method = c("AlCoV-LM", "AlCoV-Robust", "AlCoV-NNLS")) {
    if(any(coco$coverage == 0)) {
        varmat <- varmat[, coco$coverage != 0]
        coco <- coco[coco$coverage != 0,]
    }
    freqs <- coco$count/coco$coverage
    df <- as.data.frame(t(varmat))
    df$freq <- freqs

    if(method[1] == "AlCoV-LM") {
        res <- summary(lm(freq ~ 0 + ., data = df))$coef[, 1:2]
        res <- data.frame(variant = rownames(res), rho = res[,1], se = res[,2], 
            method = method)
        row.names(res) <- NULL
    } else if(method[1] == "AlCoV-Robust"){
        res <- summary(MASS::rlm(freq ~ 0 + ., data = df))$coef[, 1:2]
        res <- data.frame(variant = rownames(res), rho = res[,1], se = res[,2], 
            method = method)
        row.names(res) <- NULL
    } else if(method[1] == "AlCoV-NNLS"){
        freqs[freqs == 0] <- 0.0001
        res <- data.frame(variant = rownames(varmat),
            rho = nnls(t(varmat), freqs)$x, 
            method = method)
        row.names(res) <- NULL
    } 
    res
}

freyja <- function(coco, varmat) {
    rho_init <- provoc:::rho_initializer(varmat)

    objective <- function(rho, count, varmat, coverage) {
        freq <- count/coverage
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        prob[prob == 0 & count != 0] <- 0.000001
        if(any(is.na(freq))) {
            prob <- prob[!is.na(freq)]
            freq <- freq[!is.na(freq)]
        }
        sum(abs(freq - prob))
    }

    # Constraints will kill me --------------------------------
    # sum(p) < 1 => sum(p) - 1 < 0 => -sum(p) + 1 > 0
    u_sum1 <- rep(-1, length(rho_init))
    c_sum1 <- -1

    # p_i > 0 => 1p_i + 0p_j > 0
    u_p0 <- diag(length(rho_init))
    c_p0 <- rep(0, length(rho_init))

    ui <- rbind(u_sum1, u_p0)
    ci <- c(c_sum1, c_p0)


    res <- stats::constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        count = coco$count, coverage = coco$coverage, varmat = varmat,
        control = list(maxit = 10000))
    data.frame(variant = rownames(varmat), rho = res$par, method = "Freyja")
}

avg_freq <- function(fused, method = c("Simple Avg", "binomial", "quasibinomial")) {
    vars <- which(startsWith(colnames(fused), "var_"))
    avg <- c()
    se <- c()
    for(i in vars) {
        others <- vars[vars != i]
        sub_fuse <- fused[fused[,i] == 1 &
            apply(fused[, others, drop = FALSE], 1, sum) == 0,]
        if(method[1] == "Simple Avg") {
            if(nrow(sub_fuse) < 1) {
                avg <- c(avg, 0)
                se <- c(se, NA)
            } else {
                counts <- sub_fuse$count
                covers <- sub_fuse$coverage
                avg_tmp <- counts/covers
                avg_tmp[covers == 0] <- 0
                avg <- c(avg, mean(avg_tmp))
                se <- c(se, NA)
            }
        } else {
            res <- glm(cbind(count, coverage - count) ~ 1, 
                data = sub_fuse, 
                family = ifelse(method[1] == "binomial", "binomial", "quasibinomial"))
            avg <- c(avg, res$coef[1])
            se <- c(se, summary(res)$coef[1,2])
        }
    }
    if(sum(avg) > 1) avg <- avg / sum(avg)
    data.frame(variant = gsub("var_", "", fixed = TRUE, names(fused[, vars])),
        rho = avg,
        se = se,
        method = method)
}

provoc_optim2 <- function(coco, varmat) {
    res <- tryCatch(provoc_optim(coco, varmat)$res_df, 
        error = function(e) e)
    if("error" %in% class(res)) {
        res <- tryCatch(provoc_optim(coco, varmat)$res_df, 
            error = function(e) e)
        if("error" %in% class(res)) {
            res <- data.frame(variant = rownames(varmat), rho = NA)
        }
    }
    res$method <- "ProVoC"
    res
}

scaled_squared_counts <- function(coco, varmat) {
    rho_init <- provoc:::rho_initializer(varmat)

    objective <- function(rho, count, varmat, coverage) {
        freq <- count/coverage
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        prob[prob == 0 & count != 0] <- 0.00001
        if(any(is.na(freq))) {
            prob <- prob[!is.na(freq)]
            count <- count[!is.na(freq)]
            coverage <- coverage[!is.na(freq)]
            freq <- freq[!is.na(freq)]
        }
        preds <- prob * coverage
        denom <- coverage*prob*(1-prob)
        denom[preds-count == 0] <- 1
        sum((preds - count)^2/(denom))
    }

    # Constraints will kill me --------------------------------
    # sum(p) < 1 => sum(p) - 1 < 0 => -sum(p) + 1 > 0
    u_sum1 <- rep(-1, length(rho_init))
    c_sum1 <- -1

    # p_i > 0 => 1p_i + 0p_j > 0
    u_p0 <- diag(length(rho_init))
    c_p0 <- rep(0, length(rho_init))

    ui <- rbind(u_sum1, u_p0)
    ci <- c(c_sum1, c_p0)


    res <- stats::constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        count = coco$count, coverage = coco$coverage, varmat = varmat,
        control = list(maxit = 10000))
    data.frame(variant = rownames(varmat), rho = res$par, method = "Scaled_Squared_Counts")
}

squared_counts <- function(coco, varmat) {
    rho_init <- provoc:::rho_initializer(varmat)

    objective <- function(rho, count, varmat, coverage) {
        freq <- count/coverage
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        prob[prob == 0 & count != 0] <- 0.00001
        if(any(is.na(freq))) {
            prob <- prob[!is.na(freq)]
            count <- count[!is.na(freq)]
            coverage <- coverage[!is.na(freq)]
            freq <- freq[!is.na(freq)]
        }
        preds <- prob * coverage
        sum((preds - count)^2)
    }

    # Constraints will kill me --------------------------------
    # sum(p) < 1 => sum(p) - 1 < 0 => -sum(p) + 1 > 0
    u_sum1 <- rep(-1, length(rho_init))
    c_sum1 <- -1

    # p_i > 0 => 1p_i + 0p_j > 0
    u_p0 <- diag(length(rho_init))
    c_p0 <- rep(0, length(rho_init))

    ui <- rbind(u_sum1, u_p0)
    ci <- c(c_sum1, c_p0)


    res <- stats::constrOptim(rho_init,
        f = objective, grad = NULL,
        ui = ui, ci = ci,
        count = coco$count, coverage = coco$coverage, varmat = varmat,
        control = list(maxit = 10000))
    data.frame(variant = rownames(varmat), rho = res$par, method = "Squared_Counts")
}


truth <- list(
    Deltacron = data.frame(
        variants = c("BA.1", "BA.2", "B.1.1.529", "B.1.617.2", "B.1.617.2+K417N"),
        probs = c(0.3, 0.2, 0.1, 0.35, 0.05)
    ),
    Delta_Plus_Plus = data.frame(
        variants = c("B.1.617.2", "B.1.617.2+K417N"),
        probs = c(0.97, 0.3)
    ),
    Too_Close_to_Call = data.frame(
        variants = c("B.1.617.2", "B.1.617.2+K417N", "A.23.1", "A.23.1+E484K", 
            "BA.1", "BA.2"),
        probs = c(0.5, 0.03, 0.2, 0.02, 0.2, 0.05)
    ),
    Total_Below_One = data.frame(
        variants = c("BA.1", "BA.2", "B.1.617.2"),
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
        # Fuse to ensure same mutation list with correct order
        cocovar <- provoc:::fission(fused)
        coco <- cocovar$coco
        varmat2 <- cocovar$varmat
        all_res_tmp <- bind_rows(
            alcov(coco, varmat2, method = "AlCoV-LM"),
            tryCatch(alcov(coco, varmat2, method = "AlCoV-Robust"), 
                error = function(e) data.frame(rho = NA)),
            alcov(coco, varmat2, method = "AlCoV-NNLS"),
            freyja(coco, varmat2),
            avg_freq(fused, method = "Simple Avg"),
            #avg_freq(fused, method = "binomial"),
            #avg_freq(fused, method = "quasibinomial"),
            provoc_optim2(coco, varmat2),
            scaled_squared_counts(coco, varmat2),
            squared_counts(coco, varmat2)
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
