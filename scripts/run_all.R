suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    theme_set(theme_bw())
    library(provoc)
    library(lubridate)
    library(here)
})

animal <- read.csv(here("data/clean", "nml.csv")) 
coco <- animal %>%
    mutate(location = ifelse(grepl("MMN", sample), 
        yes = "MMN", 
        no = ifelse(grepl("MMS", sample),
            yes = "MMS",
            no = "Other")),
        mutation = label,
        coverage = alt_dp + ref_dp, # TODO: Is this correct for coverage?
        count = round(frequency * coverage, 0))

varmat <- astronomize()

fused <- fuse(coco, varmat)




#### OPTIM ------------------------------------------------
res1 <- provoc(fused, method = "optim")
any(sapply(res1, function(x) x$convergence)) # FALSE is good





#### RUNJAGS ----------------------------------------------
res2 <- provoc(fused, method = "runjags")
any(sapply(res2, function(x) x$convergence))





#### ALL TOGETHER NOW -------------------------------------
for(i in seq_along(res1)) {
    dfa <- bind_cols(res1[[i]]$point_est, res1[[i]]$sample_info)
    dfa$method <- "optim"
    dfb <- bind_cols(res2[[i]]$point_est, res2[[i]]$sample_info)
    dfb$method <- "runjags"

    if(i == 1) {
        df <- bind_rows(dfa, dfb)
    } else {
        df <- bind_rows(df, dfa, dfb)
    }
}

ggplot(df) + 
    aes(x = date, y = rho, shape = method, colour = location) +
    geom_point() + 
    geom_line(aes(group = paste(method, location)), se = FALSE) +
    facet_wrap(~ variant) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggplot(df) +
    aes(x = date, y = rho, fill = variant) +
    geom_col(colour = 1) +
    facet_wrap(~ method + location) +
    coord_flip()
