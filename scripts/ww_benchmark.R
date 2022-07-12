t0 <- Sys.time()

library(provoc)
library(dplyr)
library(ggplot2)
library(here)
library(jsonlite)

mapped <- list.files(here("output/gromstole"), 
    pattern = "mapped", full.names = TRUE)
coverage <- list.files(here("output/gromstole"), 
    pattern = "coverage", full.names = TRUE)

length(mapped)
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

t_intake <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Data intake took ", t_intake))
t1 <- Sys.time()

m_rep <- readRDS(here("output", "mutations.RDS"))
varmat <- varmat_from_variants(consensuses, max_n = 400, top_quantile = 0, mutations = m_rep, mutation_format = "aa")
varmat <- varmat[, colSums(varmat) > 1]
varmat <- varmat[rowSums(varmat) > 0, ]

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


mut_bad <- which(is.na(match(colnames(varmat), coco$mutation)))
sort(colnames(varmat)[mut_bad])


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


t_fusion1 <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Varmat creation and fusion took ", t_fusion1))
t1 <- Sys.time()

res1 <- provoc(fused = fused, method = "optim")

t_res1 <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Fitting took ", t_res1))
t1 <- Sys.time()

res1 %>% filter(rho > 0.05) %>%
    ggplot() +
        aes(x = variant, y = rho) +
        geom_point() +
        facet_wrap(~sample, scales = "free_y") +
        coord_flip()

res1 %>% 
    group_by(sample) %>%
    summarise(s = sum(rho)) %>% pull(s)

res1$mutation_list <- "Processed from NextStrain"



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


t_fusion2 <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Varmat creation and fusion took ", t_fusion2))
t1 <- Sys.time()

res2 <- provoc(fused = f2, method = "optim")
res2$mutation_list <- "Constellations repo"

t_res2 <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Fitting took ", t_res2))
t1 <- Sys.time()

ggplot(res2) +
    aes(x = variant, y = rho) +
    geom_point() +
    facet_wrap(~sample)

res2 %>% 
    group_by(sample) %>%
    summarise(s = sum(rho)) %>% pull(s)

ww_bench <- rbind(res1, res2)
ww_bench[, c("ci_low", "ci_high")] <- NULL
names(ww_bench)[1] <- "proportion"

write.csv(ww_bench, here("output", "ww_results.csv"))

print(paste0("Data intake took ", t_intake))
print(paste0("Variant 1 took ", t_fusion1))
print(paste0("Fitting 1 took ", t_res1))
print(paste0("Variant 2 took ", t_fusion2))
print(paste0("Fitting 2 took ", t_res2))
print(paste0("Total Rscript: ", difftime(Sys.time(), t0, units = "secs")))

varmat <- astronomize()
varmat <- varmat[rownames(varmat) %in% c("B.1.617.2", 
    "BA.1", "BA.2", "AY.4", "AY.4.2"),]
varmat <- varmat[, colSums(varmat) > 0]
colnames(varmat)[which(colnames(varmat) == "+2205.GAGCCAGAA")] <- "ins:22205:9"
colnames(varmat)[which(colnames(varmat) == "+8262.AACA")] <- "ins:28262:4"
colnames(varmat)[which(colnames(varmat) == "28271-")] <- "del:28271:1"
# Spike protein starts at position 21562
# Position reported as index of amino acids, hence 3*246
colnames(varmat)[which(colnames(varmat) == "aa:S:RSYLTPG246-")] <- paste0("del:", 21562 + 3*246, ":21")
colnames(varmat)[which(colnames(varmat) == "aa:S:Y144-")] <- paste0("del:", 21562 + 3*144, ":1")
colnames(varmat)[which(colnames(varmat) == "aa:S:HV69-")] <- paste0("del:", 21562 + 3*69, ":2")
# ORF 1a starts at 265
colnames(varmat)[which(colnames(varmat) == "aa:orf1a:SGF3675-")] <- paste0("del:", 265 + 3*3675 - 1, ":9")

coco <- full_join(mfiles, cfiles, by = c("position", "sample"))
coco$coverage <- apply(coco[, c("coverage.x", "coverage.y")], 1, 
    mean, na.rm = TRUE)
coco$coverage.x <- coco$coverage.y <- NULL
coco$count <- coco$frequency * coco$coverage


mut_bad <- which(is.na(match(colnames(varmat), coco$mutation)))
mut_bad <- mut_bad[!grepl(":", colnames(varmat)[mut_bad])]
for(i in mut_bad){
    col_bad <- colnames(varmat)[i]
    col_good <- parse_mutation("~", 
        substr(col_bad, 1, nchar(col_bad) - 1),
        substr(col_bad, nchar(col_bad), nchar(col_bad)))
    colnames(varmat)[i] <- col_good

}

mut_bad <- which(is.na(match(colnames(varmat), coco$mutation)))
sort(colnames(varmat)[mut_bad])

fused <- fuse(coco, varmat)

f2 <- fused %>%
    group_by(sample) %>%
    mutate(keep = n() > 15) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-keep)


t_fusion2 <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Varmat creation and fusion took ", t_fusion2))
t1 <- Sys.time()

res3 <- provoc(fused = f2, method = "optim")
res3$mutation_list <- "Delta and Omicron Constellations"

t_res3 <- difftime(Sys.time(), t1, units = "secs")
print(paste0("Fitting took ", t_res3))

res3[, c("ci_low", "ci_high")] <- NULL
head(res3, 18)

ggplot(res3, aes(x = variant, y = rho, colour = variant)) + 
    geom_point() + 
    facet_wrap(~sample) + 
    coord_flip()

write.csv(res3, here("output", "ww-delta_and_omicron.csv"))


varmat <- astronomize(path = here("data", "gromstollations"))
# Manual fix
colnames(varmat)[which(colnames(varmat) == "+2205.GAGCCAGAA")] <- "ins:22205:9"
coco <- full_join(mfiles, cfiles, by = c("position", "sample"))
coco$coverage <- apply(coco[, c("coverage.x", "coverage.y")], 1, 
    mean, na.rm = TRUE)
coco$coverage.x <- coco$coverage.y <- NULL
coco$count <- coco$frequency * coco$coverage


mut_bad <- which(is.na(match(colnames(varmat), coco$mutation)))
for(i in mut_bad){
    col_bad <- colnames(varmat)[i]
    col_good <- parse_mutation("~", 
        substr(col_bad, 1, nchar(col_bad) - 1),
        substr(col_bad, nchar(col_bad), nchar(col_bad)))
    colnames(varmat)[i] <- col_good

}

fused <- fuse(coco, varmat)

# Remove samples with too few Omicron mutations
f2 <- fused %>%
    group_by(sample) %>%
    mutate(keep = n() > 15) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-keep)


t_fusion2 <- difftime(Sys.time(), t0, units = "secs")
print(paste0("Varmat creation and fusion took ", t_fusion2))
t1 <- Sys.time()

res4 <- provoc(fused = f2, method = "optim")
res4$mutation_list <- "Delta and Omicron Gromstollations"

t_res4 <- difftime(Sys.time(), t1, units = "secs")
print(paste0("Fitting took ", t_res4))

res4[, c("ci_low", "ci_high")] <- NULL
head(res4, 18)

ggplot(res4, aes(x = variant, y = rho, colour = variant)) + 
    geom_point() + 
    facet_wrap(~sample) + 
    coord_flip()

full_join(res3, res4, by = c("variant", "sample")) %>%
    ggplot() + 
        aes(x = variant, colour = variant) +
        geom_point(aes(y = rho.x), shape = 16) +
        geom_point(aes(y = rho.y), shape = 17) +
        geom_segment(aes(xend = variant, y = rho.x, yend = rho.y)) +
        facet_wrap(~sample) +
        coord_flip()

bind_rows(res3, res4) %>%
    ggplot() +
        aes(x = variant, y = rho, colour = variant, shape = mutation_list) +
        geom_point() +
        facet_wrap(~ sample) +
        coord_flip()

write.csv(res4, file = here("output", "ww_results_gromstollations.csv"))


varmat <- astronomize()
varmat <- varmat[rownames(varmat) %in% c("B.1.617.2", 
    "BA.1", "BA.2", "AY.4", "AY.4.2"),]
varmat <- varmat[, colSums(varmat) > 0]
colnames(varmat)[which(colnames(varmat) == "+2205.GAGCCAGAA")] <- "ins:22205:9"
colnames(varmat)[which(colnames(varmat) == "+8262.AACA")] <- "ins:28262:4"
colnames(varmat)[which(colnames(varmat) == "28271-")] <- "del:28271:1"
# Spike protein starts at position 21562
# Position reported as index of amino acids, hence 3*246
colnames(varmat)[which(colnames(varmat) == "aa:S:RSYLTPG246-")] <- paste0("del:", 21562 + 3*246, ":21")
colnames(varmat)[which(colnames(varmat) == "aa:S:Y144-")] <- paste0("del:", 21562 + 3*144, ":1")
colnames(varmat)[which(colnames(varmat) == "aa:S:HV69-")] <- paste0("del:", 21562 + 3*69, ":2")
# ORF 1a starts at 265
colnames(varmat)[which(colnames(varmat) == "aa:orf1a:SGF3675-")] <- paste0("del:", 265 + 3*3675 - 1, ":9")

garmat <- astronomize(path = here("data", "gromstollations"))
# Manual fix
colnames(garmat)[which(colnames(garmat) == "+2205.GAGCCAGAA")] <- "ins:22205:9"

varmat <- varmat[!rownames(varmat) %in% rownames(garmat),]
varmat <- bind_rows(as.data.frame(varmat), as.data.frame(garmat))
varmat <- as.matrix(varmat)
varmat[is.na(varmat)] <- 0
varmat <- varmat[, colSums(varmat) > 0]
varmat[, 1:5]

mut_bad <- which(is.na(match(colnames(varmat), coco$mutation)))
mut_bad <- mut_bad[!grepl(":", colnames(varmat)[mut_bad])]
for(i in mut_bad){
    col_bad <- colnames(varmat)[i]
    col_good <- parse_mutation("~", 
        substr(col_bad, 1, nchar(col_bad) - 1),
        substr(col_bad, nchar(col_bad), nchar(col_bad)))
    colnames(varmat)[i] <- col_good

}

coco <- full_join(mfiles, cfiles, by = c("position", "sample"))
coco$coverage <- apply(coco[, c("coverage.x", "coverage.y")], 1, 
    mean, na.rm = TRUE)
coco$coverage.x <- coco$coverage.y <- NULL
coco$count <- coco$frequency * coco$coverage

fused <- fuse(coco, varmat)

# Remove samples with too few Omicron mutations
f2 <- fused %>%
    group_by(sample) %>%
    mutate(keep = n() > 15) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-keep)

res5 <- provoc(fused = f2, method = "optim")
res5$mutation_list <- "Delta and Omicron Gromstollations 2"

res5[, c("ci_low", "ci_high")] <- NULL
bind_rows(res3, res4, res5) %>%
    ggplot() +
        aes(x = variant, y = rho, colour = variant, shape = mutation_list) +
        geom_point() +
        facet_wrap(~ sample) +
        coord_flip()

write.csv(bind_rows(res3, res4, res5), file = here("output", "ww_results_gromstollations2.csv"))
