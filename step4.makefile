LD := ldblocks/EUR/fourier_ls-all.bed
CHR := $(shell seq 1 22)
K := 100
CK := 5 30

all:

twas_jobs := $(foreach kk, $(K), $(foreach cc, $(CK), $(foreach tt, Q CC, $(foreach chr, $(CHR), jobs/step4/fgwas_nn/$(tt)/$(kk)_$(cc)/$(chr)-gtex.txt.gz))))

twas: $(twas_jobs)

twas-long: $(twas_jobs:.txt.gz=.long.gz)

# % = fgwas_nn/CC/100_30/1
jobs/step4/%-gtex.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	ls -1 result/$*/*.snp-factor.gz | xargs -I file basename file .snp-factor.gz | awk -vGWAS=$(shell echo $* | awk -F'/' '{ print $$1 FS $$2 FS $$3 }') '{ print "./make.run-twas-gtex-v6.R" FS $$1 FS ("result/" GWAS) FS ("twas/gtex-v6/" GWAS "/" $$1 ".txt.gz") }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N UKBB_TWAS_$(shell echo $* | sed 's/\//_/g') -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/step4/%.long.gz: jobs/step4/%.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip > $@
	zcat $< | awk 'system("! [ -f " $$NF " ]") == 0' | gzip >> $@
	if [ $$(zcat $@ | wc -l) -gt 0 ]; then qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=36:00:00 -b y -j y -N UKBB_LONG_TWAS -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@; fi

