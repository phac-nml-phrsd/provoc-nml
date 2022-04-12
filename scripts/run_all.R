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
convergence(res1)





#### RUNJAGS ----------------------------------------------
res2 <- provoc(fused, method = "runjags")
convergence(res2)





#### ALL TOGETHER NOW -------------------------------------
res1$method <- "optim"
res2$method <- "runjags"
df <- bind_rows(res1, res2)

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
