all:

LD := ldblocks/nygcresearch-ldetect-data-ac125e47bf7f/EUR/fourier_ls-all.bed
RDIR := result/20180608
nLD := $(shell cat $(LD) | tail -n+2 | wc -l)
DATA := data/

PHENO_FILE := phenotypes/ukbb_pheno_combined.txt

queue_facto: $(foreach k, 350, jobs/20180608/facto-$(k).txt.gz)

long_facto: $(foreach k, 350, jobs/20180608/facto-$(k).long.gz)

jobs/20180608/facto-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.factorization.R" FS $$1 FS "$(LD)" FS "$(DATA)" FS "$(PHENO_FILE)" FS $* FS ("$(RDIR)/factorization/$*/" $$1) }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=3g -l h_rt=12:00:00 -b y -j y -N UKBB_FACT_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/20180608/facto-%.long.gz: jobs/20180608/facto-%.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk 'system(" ! [ -f " $$NF ".conf-lodds.gz ] ") == 0' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=6g -l h_rt=36:00:00 -b y -j y -N UKBB_long_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################
queue_conf: $(foreach k, 350, jobs/20180608/conf-$(k).txt.gz)

jobs/20180608/conf-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.confounder.R" FS $$1 FS "$(LD)" FS "$(RDIR)/factorization/$*/" FS ("$(RDIR)/confounder/$*/" $$1 ".txt.gz") }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=6g -l h_rt=4:00:00 -b y -j y -N UKBB_CONF_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@


################################################################
queue_fgwas: $(foreach k, 350, jobs/20180608/fgwas_nn-$(k).txt.gz jobs/20180608/fgwas-$(k).txt.gz)

queue_fgwas_long: $(foreach k, 350, jobs/20180608/fgwas_nn-$(k).long.gz jobs/20180608/fgwas-$(k).long.gz)

jobs/20180608/fgwas%.long.gz: jobs/20180608/fgwas%.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk 'system(" ! [ -f " $$8 ".resid.gz ] ") == 0' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=36:00:00 -b y -j y -N UKBB_long$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/20180608/fgwas-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.fgwas.R" FS $$1 FS "$(LD)" FS "data" FS "$(PHENO_FILE)" FS ("$(RDIR)/confounder/$*/" $$1 ".txt.gz") FS $* FS ("$(RDIR)/fgwas/$*/" $$1) }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=3g -l h_rt=18:00:00 -b y -j y -N UKBB_FGWAS_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/20180608/fgwas_nn-%.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	seq 1 $(nLD) | awk '{ print "./make.fgwas.R" FS $$1 FS "$(LD)" FS "data" FS "$(PHENO_FILE)" FS ("$(RDIR)/confounder/$*/" $$1 ".txt.gz") FS $* FS ("$(RDIR)/fgwas_nn/$*/" $$1) FS "TRUE" }' | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=3g -l h_rt=18:00:00 -b y -j y -N UKBB_FGWAS_NN_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################
queue_summary: $(foreach k, 350, jobs/20180608/summary-$(k).txt.gz)

jobs/20180608/summary-%.txt.gz:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	echo "./make.combine-fgwas.R $(RDIR)/fgwas/$*/ $(LD) 0.9 $(RDIR)/fgwas-$*-09" | gzip >> $@
	echo "./make.combine-fgwas.R $(RDIR)/fgwas_nn/$*/ $(LD) 0.9 $(RDIR)/fgwas_nn-$*-09" | gzip >> $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=32g -l h_rt=1:30:00 -b y -j y -N UKBB_SUM_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

