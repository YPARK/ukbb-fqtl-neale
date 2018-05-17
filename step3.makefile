
TEMPDIR := ./data
CHR := $(shell seq 1 22)
K := 100
CK := 30

all:

confounder: jobs/step3-conf.txt.gz

confounder-long: jobs/step3-long-conf.txt.gz

combine-confounder: ./make.combine-confounders.R
	[ -f $@.log ] || qsub -P compbio_lab -o $@.log -cwd -V -l h_vmem=16g -l h_rt=4:00:00 -b y -j y -N $@ ./run.sh ./$<

fgwas: jobs/step3-fgwas.txt.gz

fgwas-long: jobs/step3-long-fgwas.txt.gz

################################################################
# queue jobs to run on each LD block
jobs/step3-%.txt.gz: $(foreach chr, $(CHR), $(shell [ -d $(TEMPDIR)/$(chr)/ ] && ls -1 $(TEMPDIR)/$(chr)/ | sed 's/\///' 2> /dev/null | awk '{ print "jobs/step3/Q/$(chr)/" $$1 "-%-jobs"; print "jobs/step3/CC/$(chr)/" $$1 "-%-jobs" }'))
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N UKBB_RUN_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/step3-long-%.txt.gz: jobs/step3-%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@zcat $< | awk 'system("! [ -f " $$NF ".zscore.gz ]") == 0 && system("! [ -f " $$NF ".conf-ind.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -gt 0 ] && qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=10g -l h_rt=16:00:00 -b y -j y -N UKBB_LONG_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################
# % = $(chr)/$(ld)
jobs/step3/%-conf-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@[ -f result/conf/$(CK)/$*.conf-ind.gz ] || echo ./make.run-conf.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $(CK) phenotypes/ukbb_pheno_CC.txt result/conf/CC/$(CK)/$* >> $@
	@[ -f result/conf/$(CK)/$*.conf-ind.gz ] || echo ./make.run-conf.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $(CK) phenotypes/ukbb_pheno_Q.txt result/conf/Q/$(CK)/$* >> $@

################################################################
# % = $(Q)/$(chr)/$(ld)
# observed
jobs/step3/%-fgwas-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	@[ -f result/fgwas/$(K)/$*.zscore.gz ] || echo ./make.run-fgwas.R $(shell echo $* | awk -F'/' '{ print ("$(TEMPDIR)/" $$2 "/" $$3) OFS "1KG_EUR/chr" $$2 OFS $(K) OFS "FALSE" OFS ("result/conf/" $$1 "/$(CK)/loco_" $$2 ".txt.gz") OFS ("phenotypes/ukbb_pheno_" $$1 ".txt") OFS ("result/fgwas/$(K)/$*" ) }') >> $@
	@[ -f result/fgwas/$(K)/$*.zscore.gz ] || echo ./make.run-fgwas.R $(shell echo $* | awk -F'/' '{ print ("$(TEMPDIR)/" $$2 "/" $$3) OFS "1KG_EUR/chr" $$2 OFS $(K) OFS "TRUE" OFS ("result/conf/" $$1 "/$(CK)/loco_" $$2 ".txt.gz") OFS ("phenotypes/ukbb_pheno_" $$1 ".txt") OFS ("result/fgwas_nn/$(K)/$*" ) }') >> $@

# % = $(chr)/$(ld)
# null
## jobs/step3/%-null-jobs:
##	@for kk in $(K); do [ -f result/null-fgwas_nn/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas-null.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} TRUE result/null-fgwas_nn/$${kk}/$* >> $@; done
##	@for kk in $(K); do [ -f result/null-fgwas/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas-null.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} FALSE result/null-fgwas/$${kk}/$* >> $@; done

################################################################
# combine results
combine: $(foreach kk, $(K), $(foreach tt, Q CC, result/fgwas_nn-$(kk)-$(tt).log result/fgwas-$(kk)-$(tt).log))

CUTOFF := 2.198 # pip > 0.9

# % = fgwas-$(K)
result/%.log: ./make.combine-result.R
	qsub -P compbio_lab -o $@ -cwd -V -l h_vmem=16g -l h_rt=16:00:00 -b y -j y -N $* ./run.sh ./make.combine-result.R result/$(shell echo $* | sed 's/-/\//g')/ $(CUTOFF) result/$*
