# The half-baked algorithmic idea

library(provoc)
library(dplyr) # I tried not to use it


vdt <- function(fused) {
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
        fused1 <- fused2[new_index, ]
        fused1$frequency <- fused1$frequency - rho[winner]
        fused1$frequency[fused1$frequency < 0] <- 0
        expected_count <- rho[winner] * fused1$coverage
        fused1$count <- fused1$count - expected_count
        fused1$count[fused1$count < 0] <- 0
        fused1$coverage <- fused1$coverage - expected_count
        fused1$coverage[fused1$coverage < 0] <- 0

        fused2 <- rbind(fused2[!new_index,], fused1)
        fused2[, varnames[winner]] <- rep(0, nrow(fused2))
    }
    rho
}

vdt(fused)

bootstrap_counts <- function(df) {
    case1 <- unlist(sapply(1:nrow(df), function(x) rep(x, df$count[x])))
    df1 <- df[case1,]
    df1$case <- 1

    case0 <- unlist(sapply(1:nrow(df), 
        function(x) rep(x, df$coverage[x] - df$count[x])))
    df0 <- df[case0, ]
    df0$case <- 0

    new_df <- rbind(df1, df0)

    new_df <- new_df[sample(1:nrow(new_df), nrow(new_df), TRUE),]
    new_df %>% select(-count, -coverage) %>%
        group_by(across(c(-case))) %>%
        summarise(count = sum(case), coverage = n(), .groups = "drop")
}

vba <- function(fused, var_replicates = 20, boot_replicates = 2) {
    t0 <- Sys.time()
    varnames <- names(fused)[startsWith(names(fused), "var_")]
    
    if(length(varnames) <= 2) {
        return(data.frame(rho = NA, se = NA, time = NA, variant = varnames, method = "vba"))
    }
    reps <- replicate(var_replicates, {
        var_rm <- sample(varnames, round(0.3*length(varnames), 0), FALSE)
        fused_vars <- fused[, names(fused)[!names(fused) %in% var_rm]]
        res <- replicate(boot_replicates, {
            vdt(bootstrap_counts(fused_vars))
        })
        res_na <- matrix(NA, ncol = ncol(res), nrow = length(var_rm))
        rownames(res_na) <- var_rm
        res <- rbind(res, res_na)
        res[order(rownames(res)),]
    })
    rho <- apply(reps, 1, mean, na.rm = TRUE)
    data.frame(rho = rho,
        se = apply(reps, 1, sd, na.rm = TRUE),
        variant = gsub("var_", "", names(rho), fixed = TRUE),
        method = "vba",
        time = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
}

if(FALSE) { # Testing
    library(ggplot2)
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

