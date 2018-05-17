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

## quantitative tab
out.tab <- pheno.tab %>%
    filter(N.non.missing > 1e5, is.na(N.cases), is.na(N.controls)) %>%
        select(- starts_with('PHESANT'), - warning.for.case.control)

write_tsv(out.tab, path = 'phenotypes/ukbb_pheno_Q.txt')
write.xlsx(out.tab, 'phenotypes/ukbb_pheno_Q.xlsx', sheetName = 'quantitative')


## case-control
out.tab <- pheno.tab %>%
    filter(N.non.missing > 1e5, N.cases > 1e4, is.na(warning.for.case.control)) %>%
        select(- starts_with('PHESANT'), - warning.for.case.control)        

write_tsv(out.tab, path = 'phenotypes/ukbb_pheno_CC.txt')
write.xlsx(out.tab, 'phenotypes/ukbb_pheno_CC.xlsx', sheetName = 'case_control')

