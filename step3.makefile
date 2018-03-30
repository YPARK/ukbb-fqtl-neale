
TEMPDIR := ./data
CHR := $(shell seq 1 22)
K := 100
CK := 50

all:

confounder: jobs/step3-conf.txt.gz

confounder-long: jobs/step3-long-conf.txt.gz

combine-confounder: ./make.combine-confounders.R
	[ -f $@.log ] || qsub -P compbio_lab -o $@.log -cwd -V -l h_vmem=16g -l h_rt=1:00:00 -b y -j y -N $@ ./run.sh ./$<

fgwas: jobs/step3-fgwas.txt.gz

fgwas-long: jobs/step3-long-fgwas.txt.gz

################################################################
# queue jobs to run on each LD block
jobs/step3-%.txt.gz: $(foreach chr, $(CHR), $(shell [ -d $(TEMPDIR)/$(chr)/ ] && ls -1 $(TEMPDIR)/$(chr)/ | sed 's/\///' 2> /dev/null | awk '{ print "jobs/step3/$(chr)/" $$1 "-%-jobs" }'))
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N UKBB_RUN_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/step3-long-%.txt.gz: jobs/step3-%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@zcat $< | awk 'system("! [ -f " $$NF ".zscore.gz ]") == 0 && system("! [ -f " $$NF ".conf-ind.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -gt 0 ] && qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=10g -l h_rt=16:00:00 -b y -j y -N UKBB_LONG_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

# % = $(chr)/$(ld)
jobs/step3/%-conf-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@[ -f result/conf/$(CK)/$*.conf-ind.gz ] || echo ./make.run-conf.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $(CK) result/conf/$(CK)/$* >> $@

# % = $(chr)/$(ld)
jobs/step3/%-fgwas-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
# observed
	@for kk in $(K); do [ -f result/fgwas_nn/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} TRUE result/conf/50/$*.conf-ind.gz result/conf/50/$*.conf-lodds.gz result/fgwas_nn/$${kk}/$* >> $@; done
	@for kk in $(K); do [ -f result/fgwas/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} FALSE result/conf/50/$*.conf-ind.gz result/conf/50/$*.conf-lodds.gz result/fgwas/$${kk}/$* >> $@; done
# null
	@for kk in $(K); do [ -f result/null-fgwas_nn/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas-null.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} TRUE result/null-fgwas_nn/$${kk}/$* >> $@; done
	@for kk in $(K); do [ -f result/null-fgwas/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas-null.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} FALSE result/null-fgwas/$${kk}/$* >> $@; done

################################################################
# combine results
combine: $(foreach kk, $(K), $(foreach rr, $(R), result/fgwas_nn-$(kk)_$(rr).log result/fgwas-$(kk)_$(rr).log))

CUTOFF := 0

# % = fgwas_nn-50_5
result/%.log: ./make.combine-result.R
	qsub -P compbio_lab -o $@ -cwd -V -l h_vmem=16g -l h_rt=16:00:00 -b y -j y -N $* ./run.sh ./make.combine-result.R result/$(shell echo $* | sed 's/-/\//g')/ $(CUTOFF) result/$*
