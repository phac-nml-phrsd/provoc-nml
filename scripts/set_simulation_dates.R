library(lubridate)
library(ggplot2)
library(dplyr)

variants <- read.csv("data/covid19-epiSummary-variants.csv") %>% 
    rename(variant = Variant.Grouping, identifier = X_Identifier,
        lineage = Lineage.Grouped, percent = X.CT.Count.of.Sample..,
        week = Collection..week.)

head(variants)

variants %>%
    group_by(variant, week) %>%
    mutate(percent = sum(percent)) %>%
    ggplot() + 
        aes(x = ymd(week), y = percent, colour = variant) +
        geom_point() + geom_line() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

dates <- ymd("2021-05-30") + 7*(0:10)
variants %>%
    filter(week %in% as.character(dates)) %>%
    group_by(lineage) %>%
    mutate(max_perc = rep(max(percent), n())) %>%
    ungroup() %>%
    filter(max_perc > 0.01) %>%
    select(-max_perc) %>%
    ggplot() + 
        aes(x = ymd(week), y = percent, colour = variant, shape = lineage) +
        geom_point(size = 5) + geom_line(size = 1) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        scale_shape_manual(values = 1:20) +
        scale_x_date(breaks = dates[2 * (1:10) - 1]) +
        labs(x = "Date", y = "Percent of CT Samples") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))

variants %>%
    filter(week %in% as.character(dates)) %>%
    group_by(lineage) %>%
    mutate(max_perc = rep(max(percent), n())) %>%
    ungroup() %>%
    filter(max_perc > 0.01) %>%
    select(-max_perc) %>%
    ggplot() + 
        aes(x = ymd(week), y = percent, fill = lineage) +
        geom_area(colour = 1) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        scale_shape_manual(values = 1:20) +
        scale_x_date(breaks = dates[2 * (1:10) - 1]) +
        labs(x = "Date", y = "Percent of CT Samples") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1))
