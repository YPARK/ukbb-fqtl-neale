#!/usr/bin/env Rscript
## read the data and estimate FSR using ASHR

source('util.R')
library(dplyr)
library(readr)
library(tidyr)

trait.factors <- read_tsv('result/ukbb-fqtl-tot-traits.txt.gz')

## LFSR estimation
trait.factors <- trait.factors %>%
    mutate(pip = 1/(1+exp(-lodds))) %>%
    mutate(pos.prob = pnorm(0, mean = theta, sd = theta.se)) %>%
    mutate(neg.prob = 1 - pos.prob) %>%
    mutate(lfsr = 1 - pip * pmax(pos.prob, neg.prob))

## Calibrate overall FSR rate with different PIP cutoff
take.mean.lfsr <- function(.pip) {
    ret <- trait.factors %>% filter(pip > .pip) %>%
        summarize(lfsr = mean(lfsr))
    return(data.frame(pip.cutoff = .pip, FSR = as.numeric(ret$lfsr)))
}

n.tests <- nrow(trait.factors)
pip.cutoff <- 1/(1+exp(-c(seq(-8, 4, by = 0.5))))

fsr.tab <- do.call(rbind, lapply(pip.cutoff, take.mean.lfsr))

write_tsv(fsr.tab, path = gzfile('result/ukbb-trait-fsr.txt.gz'))

library(ggplot2)
library(scales)
library(ggrepel)

logit.trans <- trans_new(".logit", function(x) log(x) - log(1-x), function(y) 1/(1+exp(-y)))
l10.trans <- trans_new(".ln10", function(x) pmax(log10(x), -8), function(y) 10^(y))

plt <-
    gg.plot(fsr.tab, aes(x = pip.cutoff, y = FSR, label = signif(FSR, 2))) +
    geom_point() + geom_line() +
    geom_text_repel(data = fsr.tab %>% filter(FSR < .2), nudge_x = .01, direction = 'x', color = 'orange', segment.color = 'gray') +
    scale_x_continuous(breaks = signif(1/(1 + exp(-seq(-8, 4))), 2), trans = logit.trans) +
    scale_y_continuous(breaks = 10^(seq(-8, 0)), trans = l10.trans) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
          legend.title = element_blank(),
          legend.position = c(.99,.99), legend.justification = c(1,1)) +
    ylab('False sign rate') + xlab('posterior inclusion probability cutoff')

ggsave(filename = 'result/fig_fsr_pip.pdf', plot = plt, width = 4, height = 4, units = 'in')
