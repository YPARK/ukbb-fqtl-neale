## Select UKBB phenotypes with the enough number of cases

pheno.file <- 'ukbiobank_summary/phenosummary_final_11898_18597.tsv'

library(dplyr)
library(readr)
library(xlsx)
options(stringsAsFactors = FALSE)

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

################################################################
## break down into categories

cat.code.tab <- read_tsv('ukb_field.txt',
                         col_names = c('category', 'code', 'category.name'),
                         col_types = 'icc_')

library(tidyr)

pheno.cat.tab <- pheno.tab %>%
    mutate(.code = Field.code) %>%
        separate('.code', c('code', '.subcode'), sep = '_') %>%
            left_join(cat.code.tab) %>%
                filter(!is.na(category.name)) %>%
                    select(-Notes, -starts_with('PHESANT'), -warning.for.case.control, -.subcode)

small.sample <- pheno.cat.tab %>%
    filter(N.non.missing < 1e4) %>%
        select(Field.code)

unbalanced.case.control <- pheno.cat.tab %>%
    filter(!is.na(N.cases)) %>%
        filter(N.cases < 1e4 | N.controls < 1e4) %>%
            select(Field.code)

pheno.cat.tab <- pheno.cat.tab %>%
    anti_join(small.sample) %>%
        anti_join(unbalanced.case.control)

for(.cat in unique(pheno.cat.tab$category.name)) {

    .out <- pheno.cat.tab %>% filter(category.name == .cat) %>%
        select(-starts_with('category'))

    .cat.file <- .cat %>% gsub(pattern = '[- ]', replacement = '_') %>%
        (function(x) paste0('phenotypes/ukbb_', x, '.txt'))

    write_tsv(.out, path = .cat.file)
}

write_tsv(pheno.cat.tab, 'phenotypes/ukbb_all_categories.txt')
