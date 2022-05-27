# Function to generate a variant matrix

library(provoc)
library(lubridate)
library(dplyr)

N <- 1000

parse_genbank <- function(mut) {
    if(grepl("^[ATCG]", mut)) {
        c("~", as.numeric(substr(mut, 2, nchar(mut)-1)), substr(mut, nchar(mut), nchar(mut)))
    } else if(grepl("-", mut)) {
        del <- as.numeric(strsplit(mut, "-")[[1]])
        c("-", del[1], del[2]-del[1] + 1)
    } else if(grepl(":", mut)) {
        ins <- strsplit(mut, ":")[[1]]
        c("+", as.numeric(ins[1]), ins[2])
    } else {
        c(NA, NA, NA)
    }
}

sample_for_varmat <- function(lineage, count) {
    if(lineage == "Other Delta") {
        mlin <- mutations_by_lineage[grepl("AY.", mutations_by_lineage$lineage), ]
    } else {
        mlin <- mutations_by_lineage[mutations_by_lineage$lineage == lineage, ]
    }
    muts <- sample(mlin$mutation, count, replace = TRUE, prob = mlin$count/sum(mlin$count))
    muts <- t(sapply(muts, parse_genbank))
    muts <- muts[complete.cases(muts),, drop = FALSE]
    muts <- sapply(1:nrow(muts), function(x){
        provoc:::parse_mutation(muts[x,1], muts[x,2], muts[x,3])})
    unique(muts)
}

sample_varmat <- function(lineages, counts) {
    samples <- lapply(1:length(lineages), function(x) sample_for_varmat(lineages[x], counts[x]))
    names(samples) <- lineages
    provoc:::varmat_from_list(samples)
}

variants <- read.csv("data/clean/var_params.csv")

for(date in unique(variants$week)) {
    day_vars <- variants[variants$week == ymd(date),]
    for(rep in 1:10) {
        sampled_varmat <- sample_varmat(lineages = day_vars$lineage, 
            counts = round(day_vars$percent*N))
        coco <- simulate_coco(sampled_varmat, verbose = FALSE, 
            rel_counts = round(day_vars$percent*N))
        varmat <- astronomize()
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
}

head(all_res)








