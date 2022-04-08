suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(provoc)
    library(lubridate)
    library(here)
})

animal <- read.csv(here("data/clean", "nml.csv"))
s1 <- animal[animal$sample == animal$sample[1],]

varmat <- astronomize()

s1$mutation <- s1$label
s1$coverage <- s1$alt_dp + s1$ref_dp
s1$count <- round(s1$frequency*s1$coverage, 0)
coco <- s1

fused <- fuse(coco, varmat)

res1 <- provoc(fused, method = "optimize")
res2 <- provoc(fused, method = "runjags")

res2b <- melt_mcmc(res2, pivot = TRUE) %>%
    group_by(name) %>%
    summarise(par = median(value), .groups = "drop")

plot(as.numeric(res2b$par), as.numeric(res1$par))
