library(provoc)
library(dplyr)
library(ggplot2)
library(here)
library(jsonlite)

mapped <- list.files(here("output/gromstole"), 
    pattern = "mapped", full.names = TRUE)
coverage <- list.files(here("output/gromstole"), 
    pattern = "coverage", full.names = TRUE)

length(mapped) == length(coverage)

mfiles <- lapply(mapped, read.csv)
cfiles <- lapply(coverage, read.csv)
for(i in seq_along(mfiles)) mfiles[[i]]$sample <- i
for(i in seq_along(cfiles)) cfiles[[i]]$sample <- i
mfiles <- bind_rows(mfiles)
cfiles <- bind_rows(cfiles)
mfiles$position <- as.numeric(mfiles$position) + 1
cfiles$position <- as.numeric(cfiles$position) + 1

consensuses <- readLines(here("../ww_benchmark/consensus_lineages.txt")) %>%
    strsplit(split = ',|\"|"', fixed = FALSE) %>%
    sapply(`[`, 3) %>%
    unique()
consensuses <- consensuses[!is.na(consensuses)]
consensuses <- unique(c(consensuses, "BA.1", "BA.1.1"))

if(FALSE) {
    llist <- list()
    for(i in seq_along(consensuses)) {
        temp_aa <- read_json(paste0("https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage=", consensuses[i]))
        foo_aa <- data.frame(t(sapply(temp_aa[[3]], as.character)))
        llist[[consensuses[i]]] <- foo_aa
    }
    str(llist)
    varmat <- varmat_from_list(lapply(llist, function(x) x[1][x[2] > 0.8]))
    varmut <- colnames(varmat)
    varmut <- lapply(varmut, function(x) {
        c(ifelse(substr(x, nchar(x), nchar(x)) == "-", "-", "~",
            substr(x, 2, nchar(x) - 1),
            substr(x, nchar(x), nchar(x))))
    })
    varpos <- integer(length(varmut))
    for (i in seq_along(varmut)) {
        if(grepl("-", varmut[i])) {
            nums <- as.numeric(strsplit(varmut[i], split = "-")[[1]])
            varmut[i] <- paste0("-", nums[1], ".", diff(nums), collapse = "")
            varpos[i] <- nums[1]
        } else if (grepl(":", varmut[i])) {
            vals <- strsplit(varmut[i], split = ":")[[1]]
            varmut[i] <- paste0("+", vals[1], ".", vals[2])
            varpos[i] <- as.integer(vals[1])
        } else {
            varpos[i] <- as.integer(gsub("[ATCG]", "", varmut[i]))
            varmut[i] <- paste0("~", substr(varmut[i], 2, 20))
        }
    }
    colnames(varmat) <- varmut
    dim(varmat)
}

varmat <- varmat_from_variants(consensuses, max_n = 300, top_quantile = 0)
dim(varmat)
rowSums(varmat)
varmut <- colnames(varmat)
varpos <- integer(length(varmut))
for (i in seq_along(varmut)) {
    if(grepl("-", varmut[i])) {
        nums <- as.numeric(strsplit(varmut[i], split = "-")[[1]])
        varmut[i] <- paste0("-", nums[1], ".", diff(nums), collapse = "")
        varpos[i] <- nums[1]
    } else if (grepl(":", varmut[i])) {
        vals <- strsplit(varmut[i], split = ":")[[1]]
        varmut[i] <- paste0("+", vals[1], ".", vals[2])
        varpos[i] <- as.integer(vals[1])
    } else {
        varpos[i] <- as.integer(gsub("[ATCG]", "", varmut[i]))
        varmut[i] <- paste0("~", substr(varmut[i], 2, 20))
    }
}
colnames(varmat) <- varmut
dim(varmat)

coco <- full_join(
    mfiles[mfiles$position %in% varpos,], 
    cfiles[cfiles$position %in% varpos,], 
    by = c("position", "sample"))
coco$frequency[is.na(coco$frequency)] <- 0
coco$coverage <- apply(coco[, c("coverage.x", "coverage.y")], 1, 
    mean, na.rm = TRUE)
coco$coverage.x <- coco$coverage.y <- NULL
coco$count <- coco$frequency * coco$coverage

coco$aa_mutation <- coco$mutation
coco$mutation <- coco$label
coco$label <- NULL
fused <- fuse(coco, varmat)
dim(fused)

res <- provoc(fused = fused, method = "optim")

ggplot(res) +
    aes(x = variant, y = rho) +
    geom_point() +
    facet_wrap(~sample)

res %>% 
    group_by(sample) %>%
    summarise(s = sum(rho)) %>% pull(s)





varmat <- astronomize()
coco <- full_join(mfiles, cfiles, by = c("position", "sample"))
coco$coverage <- apply(coco[, c("coverage.x", "coverage.y")], 1, 
    mean, na.rm = TRUE)
coco$coverage.x <- coco$coverage.y <- NULL
coco$count <- coco$frequency * coco$coverage
fused <- fuse(coco, varmat)

f2 <- fused %>%
    group_by(sample) %>%
    mutate(keep = n() > 15) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-keep)

res2 <- provoc(fused = f2, method = "optim")

ggplot(res2) +
    aes(x = variant, y = rho) +
    geom_point() +
    facet_wrap(~sample)

res2 %>% 
    group_by(sample) %>%
    summarise(s = sum(rho)) %>% pull(s)
