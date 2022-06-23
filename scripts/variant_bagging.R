# The half-baked algorithmic idea

library(provoc)
library(dplyr) # I tried not to use it


vdt <- function(fused) {
    t0 <- Sys.time()
    varnames <- names(fused)[startsWith(names(fused), "var_")]
    if(!"frequency" %in% names(fused)) {
        fused$frequency <- fused$count / fused$coverage
    }
    fused2 <- fused
    rho <- rep(0, length(varnames))
    names(rho) <- varnames

    for(iter in seq_len(length(varnames))) {
        infos <- sapply(varnames, function(x){
            mean(fused2$frequency[fused2[, x] == 1]) - 
                mean(fused2$frequency[fused2[, x] != 1])
        })

        winner <- which.max(infos)
        if(length(winner) > 1) winner <- winner(sample(winner, 1))
        rho[winner] <- mean(fused2$frequency[fused2[, varnames[winner]] == 1])

        # Update step
        new_index <- fused2[, varnames[winner]] == 1
        expected_count <- new_index * rho[winner] * fused2$coverage
        fused2$frequency <- fused2$frequency - new_index * rho[winner]
        fused2$count <- fused2$count - expected_count
        fused2$coverage <- fused2$coverage - expected_count

        fused2$frequency[fused2$frequency < 0] <- 0
        fused2$count[fused2$count < 0] <- 0
        fused2$coverage[fused2$coverage < 0] <- 0

        fused2[, varnames[winner]] <- rep(0, nrow(fused2))
    }
    res_df <- data.frame(rho = rho,
        se = NA,
        variant = gsub("var_", "", names(rho), fixed = TRUE),
        method = "vdt",
        time = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    rownames(res_df) <- NULL
    res_df
}

bootstrap_poisson <- function(df) {
    if("frequency" %in% names(df)) {
        df$coverage <- rpois(nrow(df), df$coverage)
        df$count <- rbinom(nrow(df), prob = df$frequency, size = df$count)
        df$frequency <- df$count/df$coverage
        df$frequency[df$coverage == 0] <- 0
        return(df)
    } else {
        df$coverage <- rpois(nrow(df), df$coverage)
        df$count <- rbinom(nrow(df), prob = df$count/df$coverage, size = df$count)
        df$count[df$coverage == 0] <- 0
        return(df)
    }
}

vba <- function(fused, var_replicates = 20, boot_replicates = 2) {
    t0 <- Sys.time()
    varnames <- names(fused)[startsWith(names(fused), "var_")]
    
    if(length(varnames) <= 2) {
        return(data.frame(rho = NA, se = NA, time = NA, 
            variant = gsub("var_", "", varnames, fixed = TRUE), 
            method = "vba"))
    }
    reps <- do.call(rbind, replicate(var_replicates, {
        var_rm <- sample(varnames, round(0.3*length(varnames), 0), FALSE)
        fused_vars <- fused[, names(fused)[!names(fused) %in% var_rm]]
        res <- do.call(rbind, replicate(boot_replicates,
            vdt(bootstrap_poisson(fused_vars)), simplify = FALSE))
    }, simplify = FALSE))
    res_df <- do.call(rbind, lapply(split(reps, reps$variant), function(x) {
        data.frame(rho = mean(x$rho), se = sd(x$rho), variant = x$variant[1],
            method = "vba", time = sum(x$time))
    }))
}

if(FALSE) { # Testing
    library(ggplot2)
    library(microbenchmark)
    varmat <- astronomize()
    true_vars <- c("B.1.1.529", "BA.1", "BA.2", "B.1.617.2", "AY.4")
    true_counts <- c(500, 400, 100, 500, 100)
    varmat2 <- varmat[rownames(varmat) %in% true_vars,]
    coco <- simulate_coco(varmat2, rel_counts = true_counts)
    truth <- as.data.frame(t(matrix(true_counts / sum(true_counts))))
    names(truth) <- paste0("var_", true_vars)

    fused <- fuse(coco, varmat)
    fused <- fused[fused$coverage > 0, ]

    foo <- vba(fused, 5, 5)

    bind_rows(
        as.data.frame(t(apply(foo, c(1), mean, na.rm = TRUE))),
        as.data.frame(t(vdt(fused))),
        truth
    ) %>%
        mutate(method = c("vba", "vdt", "truth")) %>%
        tidyr::pivot_longer(-method) %>%
        ggplot() + 
            aes(x = name, y = value, colour = method) +
            geom_point() +
            coord_flip()
}

