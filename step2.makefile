################################################################
# Select UKBB phenotypes with the enough number of cases
# from the Ben Neale's data
LD := ldblocks/EUR/fourier_ls-all.bed
TEMPDIR := ./data
UKBB_Neale := ukbiobank_summary

all: phenotypes/ukbb_pheno_CC.txt phenotypes/ukbb_pheno_Q.txt phenotypes/ukbb_pheno_combined.txt

phenotypes/ukbb_pheno_Q.txt: phenotypes/ukbb_pheno_CC.txt
	echo $@

phenotypes/ukbb_pheno_CC.txt:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	R --vanilla < make.ukbb-pheno.R

phenotypes/ukbb_pheno_combined.txt: phenotypes/ukbb_pheno_CC.txt phenotypes/ukbb_pheno_Q.txt
	cat $^ | awk 'NR == 1 || $$1 != "Field.code"' > $@

PHENONAMES := $(shell [ -f phenotypes/ukbb_pheno.txt ] && tail -n+2 phenotypes/ukbb_pheno.txt | cut -f1)

distribute: jobs/step2.txt.gz

jobs/step2.txt.gz: $(foreach pheno, $(PHENONAMES), jobs/step2/$(pheno)_job.txt.gz)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=16g -l h_rt=2:00:00 -b y -j y -N UKBB_DATA -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

# distribute the data: % = 1618
jobs/step2/%_job.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@[ $(shell echo $* | wc -w) -gt 0 ] && awk 'BEGIN { print "./make.data-ldblock.R" FS "$(UKBB_Neale)/$*.assoc.tsv.gz" FS "$(TEMPDIR)/" }' | gzip > $@
