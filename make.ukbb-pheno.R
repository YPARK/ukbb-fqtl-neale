## Select UKBB phenotypes with the enough number of cases

pheno.file <- 'ukbiobank_summary/phenosummary_final_11898_18597.tsv'

library(dplyr)
library(readr)
library(xlsx)

.codes <- list.files(path = 'ukbiobank_summary/', pattern = 'assoc.tsv.gz') %>%
    sapply(gsub, pattern = '.assoc.tsv.gz', replacement = '') %>%
        unlist(use.names = FALSE)

pheno.tab <- data.frame(Field.code = .codes) %>%
    left_join(read_tsv(pheno.file))

quant.tab <- pheno.tab %>% filter(N.non.missing > 1e5,
                                  is.na(N.cases), is.na(N.controls))

case.control.tab <- pheno.tab %>% filter(N.non.missing > 1e5, N.cases > 5e4,
                                         is.na(warning.for.case.control))

out.tab <- bind_rows(quant.tab, case.control.tab) %>%
    select(- starts_with('PHESANT'), - warning.for.case.control)

write_tsv(out.tab, path = 'phenotypes/ukbb_pheno.txt')
write.xlsx(out.tab, 'phenotypes/ukbb_pheno.xlsx', sheetName = 'phenotypes')

