# Functions for fitting other peoples' models



alcov <- function(coco, varmat, method = c("lm", "robust", "binomial")) {
    cocovar <- provoc:::fission(fused)
    coco <- cocovar$coco
    varmat <- cocovar$varmat
    freqs <- coco$count/coco$coverage
    df <- as.data.frame(t(provoc:::fission(fused)$varmat))
    df$freq <- freqs
    if(method[1] == "lm") {
        res <- summary(lm(freq ~ 0 + ., data = df))$coef[, 1:2]
    } else if(method[1] == "robust"){
        res <- summary(MASS::rlm(freq ~ 0 + ., data = df))$coef[, 1:2]
    } else {
        df$freq <- NULL
        df$count <- coco$count
        df$coverage <- coco$coverage
        res <- summary(glm(cbind(count, coverage - count) ~ 0 + ., data = fused[, -3], family = binomial))$coef
    }
    res
}

freyja <- function(coco, varmat) {
    rho_init <- provoc:::rho_initializer(varmat)

    objective <- function(rho, count, varmat, coverage) {
        prob <- as.numeric(rho %*% varmat)
        prob[coverage == 0] <- 0
        prob[prob == 0 & count != 0] <- 0.000001
        sum(abs(count/coverage - prob))
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
    data.frame(variant = rownames(varmat), rho = res$par)
}

avg_freq <- function(fused) {
    vars <- which(startsWith(colnames(fused), "var_"))
    avg <- c()
    for(i in vars) {
        others <- vars[vars != i]
        sub_fuse <- fused[fused[,i] == 1 &
            apply(fused[, others], 1, sum) == 0,]
        avg <- c(avg, mean(sub_fuse$count/sub_fuse$coverage))
    }
    if(sum(avg) > 1) avg <- avg / sum(avg)
    data.frame(variant = gsub("var_", "", fixed = TRUE, names(fused[, vars])),
        rho = avg)
}


varmat <- astronomize()
varmat <- varmat[rownames(varmat) %in% c("BA.1", "BA.2", "B.1.1.529"), ]
varmat <- varmat[, apply(varmat, 2, sum) > 0]
coco <- simulate_coco(varmat)
fused <- fuse(coco, varmat)
# Fuse to ensure same mutation list with correct order
cocovar <- provoc:::fission(fused)
coco <- cocovar$coco
varmat <- cocovar$varmat

alcov(coco, varmat)
alcov(coco, varmat, method = "robust")
alcov(coco, varmat, method = "binomial")
freyja(coco, varmat)
avg_freq(fused)
provoc_optim(coco, varmat)$res_df

