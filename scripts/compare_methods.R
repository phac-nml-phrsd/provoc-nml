# Functions for fitting other peoples' models
library(nnls)


alcov <- function(coco, varmat, method = c("lm", "robust", "nnls")) {
    if(any(coco$coverage == 0)) {
        varmat <- varmat[, coco$coverage != 0]
        coco <- coco[coco$coverage != 0,]
    }
    freqs <- coco$count/coco$coverage
    df <- as.data.frame(t(varmat))
    df$freq <- freqs

    if(method[1] == "lm") {
        res <- summary(lm(freq ~ 0 + ., data = df))$coef[, 1:2]
        res <- data.frame(variant = rownames(res), rho = res[,1], se = res[,2], 
            method = method)
        row.names(res) <- NULL
    } else if(method[1] == "robust"){
        res <- summary(MASS::rlm(freq ~ 0 + ., data = df))$coef[, 1:2]
        res <- data.frame(variant = rownames(res), rho = res[,1], se = res[,2], 
            method = method)
        row.names(res) <- NULL
    } else if(method[1] == "nnls"){
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
    data.frame(variant = rownames(varmat), rho = res$par, method = "freyja")
}

avg_freq <- function(fused, method = c("simple", "binomial", "quasibinomial")) {
    vars <- which(startsWith(colnames(fused), "var_"))
    avg <- c()
    se <- c()
    for(i in vars) {
        others <- vars[vars != i]
        sub_fuse <- fused[fused[,i] == 1 &
            apply(fused[, others], 1, sum) == 0,]
        if(method[1] == "simple") {
            avg <- c(avg, mean(sub_fuse$count/sub_fuse$coverage))
            se <- c(se, NA)
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
    res$method <- "provoc"
    res
}


varmat <- astronomize()
varmat <- varmat[rownames(varmat) %in% c("BA.1", "BA.2", "B.1.1.529"), ]
varmat <- varmat[, apply(varmat, 2, sum) > 0]
true_vals <- data.frame(
    variant = row.names(varmat),
    prob = c(0.6, 0.3, 0.1)
)
rel_counts <- round(true_vals$prob * 500)

for(i in 1:1000) {
    coco <- simulate_coco(varmat, rel_counts = rel_counts, verbose = FALSE)
    fused <- fuse(coco, varmat, verbose = FALSE)
    fused <- fused[fused$coverage > 0,]
    # Fuse to ensure same mutation list with correct order
    cocovar <- provoc:::fission(fused)
    coco <- cocovar$coco
    varmat <- cocovar$varmat
    all_res_tmp <- bind_rows(
        alcov(coco, varmat, method = "lm"),
        alcov(coco, varmat, method = "robust"),
        alcov(coco, varmat, method = "nnls"),
        freyja(coco, varmat),
        avg_freq(fused, method = "simple"),
        #avg_freq(fused, method = "binomial"),
        #avg_freq(fused, method = "quasibinomial"),
        provoc_optim2(coco, varmat)
    )
    all_res_tmp$iter <- i
    if(i == 1) {
        all_res <- all_res_tmp
    } else {
        all_res <- bind_rows(all_res, all_res_tmp)
    }
}

all_res <- left_join(all_res, true_vals, by = "variant")

ggplot(all_res, aes(x = method, y = rho, fill = method)) +
    geom_violin(draw_quantiles = c(0.045, 0.5, 0.954)) +
    facet_wrap(~ variant, scales = "free_y") +
    geom_hline(aes(yintercept = prob, group = variant))

ggplot(all_res) +
    aes(x = method, y = rho - prob, fill = variant) +
    geom_violin(draw_quantiles = c(0.045, 0.5, 0.954)) +
    labs(x = NULL, y = "Error") +
    geom_hline(yintercept = 0)
