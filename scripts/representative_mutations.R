library(provoc)
library(dplyr)
library(lubridate)
library(here)

# Calculate in-group (70% of seqs in this lineage have this mutation)
m_in <- mutations_by_lineage %>%
    group_by(lineage) %>%
    mutate(q = ecdf(count)(count), max_count = max(count)) %>%
    mutate(keep = q > 0.7 & max_count > 20, in_group = count/max_count) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-max_count)

# Assumption: Each mutation appears in each lineage once
# Out-group mutations are the ones common in other lineages
mtab <- as.data.frame(table(m_in$mutation))
names(mtab) <- c("mutation", "out_group")

m_in <- left_join(m_in, mtab, by = "mutation")

# Representative - >70% in-group and <30% out-group 
m_rep <- filter(m_in, out_group < 0.3 * max(out_group))

saveRDS(m_rep, here("output", "mutations.RDS"))
