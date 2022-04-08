suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(provoc)
    library(lubridate)
    library(here)
})

animal <- read.csv(here("data/clean", "nml.csv"))
#s1 <- animal[animal$sample == animal$sample[1],]

varmat <- astronomize()

animal$mutation <- animal$label
animal$coverage <- animal$alt_dp + animal$ref_dp # TODO: Is this the correct coverage???
animal$count <- round(animal$frequency*animal$coverage, 0)
coco <- animal

fused <- fuse(coco, varmat)

res1 <- provoc(fused, method = "optim")
any(sapply(res1, function(x) x$convergence)) # FALSE is good

res1a <- lapply(seq_along(res1), function(x) {
    est <- res1[[x]]$point_est
    est$convergence <- res1[[x]]$convergence
    est$sample <- names(res1)[x]
    est
}) %>% bind_rows() 

sample_info <- animal %>%
    group_by(sample, date, region) %>%
    tally()

res1a <- left_join(res1a, sample_info, by = "sample")

ggplot(res1a, aes(x = date, y = rho)) + 
    geom_point() +
    facet_wrap(~ variant) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))




res2 <- provoc(fused, method = "runjags")
any(sapply(res2, function(x) x$convergence))

res2a <- lapply(seq_along(res2), function(x) {
    est <- res2[[x]]$point_est
    est$convergence <- res2[[x]]$convergence
    est$sample <- names(res2)[x]
    est
}) %>% bind_rows() 

res2a <- left_join(res2a, sample_info, by = "sample")

ggplot(res2a, aes(x = date, y = rho)) + 
    geom_point() +
    facet_wrap(~ variant)
