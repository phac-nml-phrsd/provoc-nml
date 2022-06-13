# Simulations for presentation

suppressPackageStartupMessages({
    library(ggplot2)
    theme_set(theme_bw())
    library(dplyr)
    library(provoc)
    library(jsonlite)
    library(here)
})


source(here("scripts", "alt_methods.R"))

format_variant <- function(variant) {
    json_file <- read_json(paste0(
        "https://lapis.cov-spectrum.org/open/v1/sample/nuc-mutations?pangoLineage=", 
        variant))
    if(length(json_file$errors) > 0) {
        print(json_file$errors)
        stop("Errors")
    }
    df <- as.data.frame(matrix(unlist(json_file$data), ncol = 3, byrow = TRUE))
    df$variant <- variant
    names(df) <- c("mutation", "frequency", "count", "variant")

    df$position <- sapply(strsplit(df$mutation, split = ""), function(x) {
        as.numeric(paste0(gsub("[A-za-z\\-]", "", x), collapse = ""))
    })
    df$alt <- sapply(strsplit(df$mutation, split = "[0-9]"),
        function(x) x[length(x)])
    df$type <- case_when(df$alt == "-" ~ "-",
        nchar(df$alt) > 1 ~ "+", TRUE ~ "~")
    df
}

sample_sequences <- function(df, N) {
    table(replicate(N, 
        df$mutation[sample(1:length(df$mutation), 
            prob = df$frequency)]))
}


if(!file.exists(here("data", "nsgb_select.csv"))) {
    variants <- c("B.1.1.529", "BA.1", "BA.2", 
        "B.1.617.2", "AY.4", "AY.4.2", "AY.25.1", "AY.103", "AY.27", "AY.74",
        "P.1", "P.1.1", "P.1.14",
        "B.1.1.7")

    df <- bind_rows(lapply(varmat, format_variant))

    write.csv(df, here("data", "nsgb_select.csv"), row.names = FALSE)
} else {
    df <- read.csv(here("data", "nsgb_select.csv"))
}



varmat <- astronomize()



