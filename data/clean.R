suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(provoc)
    library(lubridate)
    library(here)
})

nml1 <- read.csv(here("data", "raw", "VCFDATA_MontrealAug15-Oct23.csv"))
nml2 <- read.csv(here("data", "raw", "VCFDATA_VLIAug15-Oct23.csv"))

nml <- bind_rows(nml1, nml1) %>% 
    distinct()

animal <- nml %>% 
    select(
        sample = Sample, region = REGION, date = Date,
        pos = POS, alt = ALT, ref = REF, frequency = ALT_FREQ, 
        alt_dp = ALT_DP, ref_dp = REF_DP, 
        ref_aa = REF_AA, alt_codon = ALT_CODON, alt_aa = ALT_AA) %>%
    mutate(
        date = ymd(date),
        type = ifelse(grepl("\\+", alt), 
            yes = "+", 
            no = ifelse(grepl("-", alt), 
                yes = "-", 
                no = "~"))) %>%
    mutate(
        alt = gsub("-", "", gsub("+", "", alt, fixed = TRUE), fixed = TRUE))

animal$label <- sapply(seq_len(nrow(animal)), function(i) {
    parse_mutation(animal$type[i], animal$pos[i], animal$alt[i])
})

write.csv(animal, here("data", "clean", "nml.csv"), row.names = FALSE)
