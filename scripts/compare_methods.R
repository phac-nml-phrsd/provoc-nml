# Functions for fitting other peoples' models
library(nnls)
library(provoc)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

alcov <- function(coco, varmat, method = c("AlCoV-LM", "AlCoV-Robust", "AlCoV-NNLS", "AlCoV-Binom")) {
    freqs <- coco$count/coco$coverage
    df <- as.data.frame(t(varmat))
    df$freq <- freqs

    bind_rows(lapply(method, function(methi) {
        t0 <- Sys.time()
        if(methi == "AlCoV-LM") {
            res <- summary(lm(freq ~ 0 + ., data = df))$coef[, 1:2]
        } else if(methi == "AlCoV-Robust"){
            res <- tryCatch(summary(MASS::rlm(freq ~ 0 + ., data = df))$coef[, 1:2],
                error = function(e) data.frame(NA, NA))
        } else if(methi == "AlCoV-NNLS"){
            freqs[freqs == 0] <- 0.0001
            res <- cbind(nnls(t(varmat), freqs)$x, NA)
            rownames(res) <- rownames(varmat)
        } else if(methi == "AlCoV-Binom"){
            res <- summary(glm(freq ~ ., data = df))$coef[-1, 1:2]
        }
        res <- data.frame(variant = rownames(res), rho = res[,1], se = res[,2], 
            method = methi, time = as.numeric(difftime(Sys.time(), t0, units = "secs"))/nrow(res))
        row.names(res) <- NULL
        res
    }))
}

optim_methods <- function(coco, varmat, 
        method = c("freyja", "provoc", "squared", "squared_scaled")) {
    rho_init <- provoc:::rho_initializer(varmat)

    objective <- function(rho, count, varmat, coverage, methi) {
        freq <- count/coverage
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        prob[prob == 0 & count != 0] <- 0.000001
        if(any(is.na(freq))) {
            prob <- prob[!is.na(freq)]
            count <- count[!is.na(freq)]
            coverage <- coverage[!is.na(freq)]
            freq <- freq[!is.na(freq)]
        }
        if(methi == "freyja"){
            return(sum(abs(freq - prob)))
        } else if(methi == "provoc") {
            return(-sum(stats::dbinom(x = count, size = coverage, 
                prob = prob, log = TRUE)))
        } else if(methi %in% c("squared", "squared_scaled")) {
            preds <- prob * coverage
            denom <- if(methi == "squared") {
                rep(1, length(preds))
            } else { 
                coverage*prob*(1-prob)
            }
            denom[preds-count == 0] <- 1
            return(sum((preds - count)^2/(denom)))
        }
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


    bind_rows(lapply(method, function(meth) {
        t0 <- Sys.time()
        res <- stats::constrOptim(rho_init,
            f = objective, grad = NULL,
            ui = ui, ci = ci,
            count = coco$count, coverage = coco$coverage, varmat = varmat,
            methi = meth,
            control = list(maxit = 10000))
        data.frame(variant = rownames(varmat), rho = res$par, method = meth,
            time = as.numeric(difftime(Sys.time(), t0, units = "secs"))/length(res$par))
    }))
}

avg_freq <- function(fused, method = c("Simple Avg", "Simple Med", "Binomial_GLM", "Binomial_GLM_Quasi")) {
    vars <- which(startsWith(colnames(fused), "var_"))
    avg <- c()
    se <- c()
    times <- c()
    for(i in vars) {
        t0 <- Sys.time()
        others <- vars[vars != i]
        sub_fuse <- fused[fused[,i] == 1 &
            apply(fused[, others, drop = FALSE], 1, sum) == 0,]
        if(method[1] %in% c("Simple Avg", "Simple Med")) {
            if(nrow(sub_fuse) < 1) {
                avg <- c(avg, NA)
                se <- c(se, NA)
            } else {
                avg_tmp <- sub_fuse$count/sub_fuse$coverage
                avg_tmp[sub_fuse$coverage == 0] <- 0
                if(method[1] == "Simple Avg") {
                    avg <- c(avg, mean(avg_tmp))
                    se <- c(se, sd(avg_tmp))
                } else {
                    avg <- c(avg, median(avg_tmp))
                    se <- c(se, IQR(avg_tmp))
                }
            }
        } else {
            successes <- sub_fuse$count
            failures <- sub_fuse$coverage - successes
            if(length(successes) > 2) {
                res <- glm(cbind(successes, failures) ~ 1,
                    family = ifelse(method[1] == "Binomial_GLM", binomial, quasibinomial))
                avg <- c(avg, exp(res$coef)/(1 + exp(res$coef)))
                se_tmp <- summary(res)$coef[1,2]
                se <- c(se, exp(se_tmp)/(1 + exp(se_tmp)))
            } else {
                avg <- c(avg, NA)
                se <- c(se, NA)
            }
        }
        times <- c(times, as.numeric(difftime(Sys.time(), t0, units = "secs")))
    }
    #if(sum(avg) > 1) avg <- avg / sum(avg)
    data.frame(variant = gsub("var_", "", fixed = TRUE, names(fused[, vars])),
        rho = avg,
        se = se,
        method = method,
        time = times)
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
