# Process the NSGB Subsets

library(here)
library(dplyr)
library(provoc)

vectors_to_matrix <- function(list_of_mutations){
    # Only take mutations present in at least 5% of sequences
    mut_table <- table(unlist(list_of_mutations))
    all_vals <- names(mut_table[mut_table/length(list_of_mutations) > 0.05])
    out_mat <- matrix(0, ncol = length(all_vals), nrow = length(list_of_mutations))
    colnames(out_mat) <- all_vals
    for(i in seq_along(list_of_mutations)) {
        out_mat[i, list_of_mutations[[i]][list_of_mutations[[i]] %in% all_vals]] <- 1
    }
    out_mat
}

lineages <- list.files(here("output"), pattern = "samples.txt", full.names = TRUE)

for(i in seq_along(lineages)) {
    print(lineages[i])
    li <- read.csv(lineages[i], sep = ";", strip.white = TRUE)
    print(dim(li))

    sub <- vectors_to_matrix(strsplit(li$substitutions, split = ","))
    for(ii in seq_along(colnames(sub))) {
        sub_name <- colnames(sub)[ii]
        colnames(sub)[ii] <- parse_mutation(type = "~",
            pos = substr(sub_name, 2, nchar(sub_name) - 1),
            alt = substr(sub_name, nchar(sub_name), nchar(sub_name)))
    }
    ins <- vectors_to_matrix(strsplit(li$insertions, split = ","))
    if(ncol(ins) > 0) {
        for(ii in seq_along(colnames(ins))) {
            ins_name <- colnames(ins)[ii]
            inss2 <- strsplit(ins_name, ":")[[1]]
            colnames(ins)[ii] <- paste0("ins:", inss2[1], ":", nchar(inss2[2]))
        }
    }
    del <- vectors_to_matrix(strsplit(li$deletions, split = ","))
    for(ii in seq_along(colnames(del))) {
        del_name <- colnames(del)[ii]
        dels2 <- as.numeric(strsplit(del_name, "-")[[1]])
        colnames(del)[ii] <- paste0("del:", dels2[1], ":", diff(dels2))
    }

    df <- cbind(li[, c("date", "region", "country", "length", "pango_lineage", "date_submitted")], sub) %>% 
        cbind(ins) %>% 
        cbind(del)

    lin_file <- gsub("_samples.txt", "_samplemat.csv", lineages[i])
    write.csv(df, file = lin_file, row.names = FALSE)
}

