
alcov <- function(coco, varmat, method = c("AlCoV-LM", "AlCoV-Robust", "AlCoV-NNLS", "AlCoV-Binom")) {
    freqs <- coco$count/coco$coverage
    df <- as.data.frame(t(varmat))
    df$freq <- freqs
    df <- df[!is.na(df$freq),]

    bind_rows(lapply(method, function(methi) {
        t0 <- Sys.time()
        if(methi == "AlCoV-LM") {
            res <- summary(lm(freq ~ 0 + ., data = df))$coef[, 1:2, drop = FALSE]
        } else if(methi == "AlCoV-Robust"){
            res <- tryCatch(summary(MASS::rlm(freq ~ 0 + ., data = df))$coef[, 1:2, drop = FALSE],
                error = function(e) data.frame(NA, NA))
        } else if(methi == "AlCoV-NNLS"){
            df$freq[df$freq == 0] <- 0.0001
            res <- cbind(nnls(t(varmat), df$freq)$x, NA)
            rownames(res) <- rownames(varmat)
        } else if(methi == "AlCoV-Binom"){
            res <- summary(glm(freq ~ ., data = df))$coef[-1, 1:2, drop = FALSE]
        }
        res <- data.frame(variant = as.character(rownames(res)), 
            rho = as.numeric(res[,1]), 
            se = as.numeric(res[,2]), 
            method = as.character(methi), 
            time = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))/nrow(res), 4))
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
        data.frame(variant = as.character(rownames(varmat)),
            rho = as.numeric(res$par), 
            se = NA, 
            method = as.character(meth),
            time = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))/length(res$par), 5))
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
                    avg <- c(avg, median(avg_tmp, na.rm = TRUE))
                    se <- c(se, IQR(avg_tmp, na.rm = TRUE))
                }
            }
        } else {
            successes <- sub_fuse$count
            failures <- sub_fuse$coverage - successes
            if(length(successes) > 3) {
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
