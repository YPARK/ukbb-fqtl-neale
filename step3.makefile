
TEMPDIR := ./data
CHR := $(shell seq 1 22)
K := 100
CK := 5 30

all:

confounder: jobs/step3-conf.txt.gz

confounder-long: jobs/step3-long-conf.txt.gz

combine-confounder: jobs/step3/conf-combine.txt.gz

fgwas: jobs/step3-fgwas_nn.txt.gz jobs/step3-fgwas.txt.gz

fgwas-long: jobs/step3-long-fgwas.txt.gz

################################################################
## confounder jobs
jobs/step3-conf.txt.gz: $(foreach F, $(CK), $(foreach chr, $(CHR), $(shell [ -d $(TEMPDIR)/$(chr)/ ] && ls -1 $(TEMPDIR)/$(chr)/ | sed 's/\///' 2> /dev/null | awk '{ print "jobs/step3/$(chr)/" $$1 "/" $(F) "-conf-jobs" }')))
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N UKBB_CONF_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

################################################################
# % = $(chr)/$(ld)/$(cf)
# argv[1] = in.dir; argv[2] = plink.hdr; argv[3] = factors; argv[4] = phenotype; argv[4] = output
jobs/step3/%-conf-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@printf "" > $@
	echo ./make.run-conf.R $(TEMPDIR)/$(shell echo $* | awk -F'/' '{ print $$1 FS $$2 }') 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $(shell echo $* | awk -F'/' '{ print $$3 }') phenotypes/ukbb_pheno_CC.txt result/conf/CC/$(shell echo $* | awk -F'/' '{ print $$3 }')/$(shell echo $* | awk -F'/' '{ print $$1 FS $$2 }') >> $@
	echo ./make.run-conf.R $(TEMPDIR)/$(shell echo $* | awk -F'/' '{ print $$1 FS $$2 }') 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $(shell echo $* | awk -F'/' '{ print $$3 }') phenotypes/ukbb_pheno_Q.txt result/conf/Q/$(shell echo $* | awk -F'/' '{ print $$3 }')/$(shell echo $* | awk -F'/' '{ print $$1 FS $$2 }') >> $@

jobs/step3/conf-combine.txt.gz: $(foreach F, $(CK), $(foreach tt, Q CC, jobs/step3/temp-conf-combine-$(tt)-$(F)))
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $^ | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=16g -l h_rt=4:00:00 -b y -j y -N UKBB_CONF_COMBINE -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@
	rm $^

# % = $(tt)-$(F)
jobs/step3/temp-conf-combine-%:
	echo ./make.combine-confounders.R result/conf/$(shell echo $* | sed 's/-/\//') > $@

################################################################
# queue jobs to run on each LD block
jobs/step3-%.txt.gz: $(foreach chr, $(CHR), $(foreach tt, Q CC, $(foreach kk, $(K), $(foreach cc, $(CK), jobs/step3/$(chr)/$(tt)/$(kk)/$(cc)-%-jobs))))
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@cat $^ | gzip > $@
	qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=8g -l h_rt=4:00:00 -b y -j y -N UKBB_RUN_$* -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/step3-fgwas%.long.gz: jobs/step3-fgwas%.txt.gz
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@zcat $< | awk 'system("! [ -f " $$NF ".zscore.gz ]") == 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -gt 0 ] && qsub -P compbio_lab -o /dev/null -binding "linear:1" -cwd -V -l h_vmem=10g -l h_rt=16:00:00 -b y -j y -N UKBB_LONG_FGWAS -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

# observed
# % = $(chr)/$(tt)/$(kk)/$(cc)
jobs/step3/%-fgwas_nn-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	ls -1 $(TEMPDIR)/$(shell echo $* | awk -F'/' '{ print $$1 }')/ | sed 's/\///g' | awk -vchr=$(shell echo $* | awk -F'/' '{ print $$1 }') -vtt=$(shell echo $* | awk -F'/' '{ print $$2 }') -vkk=$(shell echo $* | awk -F'/' '{ print $$3 }') -vcc=$(shell echo $* | awk -F'/' '{ print $$4 }') '{ ld = $$1; print "./make.run-fgwas.R" OFS ("$(TEMPDIR)/" chr "/" ld) OFS ("1KG_EUR/chr" chr) OFS (kk) OFS "TRUE" OFS ("result/conf/" tt "/" cc "/loco_" chr ".txt.gz") OFS ("phenotypes/ukbb_pheno_" tt ".txt") OFS ("result/fgwas_nn/" tt "/" kk "_" cc "/" chr "/" ld) }' > $@

# % = $(chr)/$(tt)/$(kk)/$(cc)
jobs/step3/%-fgwas-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	ls -1 $(TEMPDIR)/$(shell echo $* | awk -F'/' '{ print $$1 }')/ | sed 's/\///g' | awk -vchr=$(shell echo $* | awk -F'/' '{ print $$1 }') -vtt=$(shell echo $* | awk -F'/' '{ print $$2 }') -vkk=$(shell echo $* | awk -F'/' '{ print $$3 }') -vcc=$(shell echo $* | awk -F'/' '{ print $$4 }') '{ ld = $$1; print "./make.run-fgwas.R" OFS ("$(TEMPDIR)/" chr "/" ld) OFS ("1KG_EUR/chr" chr) OFS (kk) OFS "FALSE" OFS ("result/conf/" tt "/" cc "/loco_" chr ".txt.gz") OFS ("phenotypes/ukbb_pheno_" tt ".txt") OFS ("result/fgwas/" tt "/" kk "_" cc "/" chr "/" ld) }' > $@


################################################################
# combine results
combine: $(foreach kk, $(K), $(foreach cc, $(CK), $(foreach tt, Q CC, result/fgwas_nn-$(tt)-$(kk)_$(cc).log result/fgwas-$(tt)-$(kk)_$(cc).log)))

CUTOFF := 2.198 # pip > 0.9

# % = fgwas-$(T)-$(K)_$(CK)
result/%.log: ./make.combine-result.R
	qsub -P compbio_lab -o $@ -cwd -V -l h_vmem=16g -l h_rt=16:00:00 -b y -j y -N $* ./run.sh ./make.combine-result.R result/$(shell echo $* | sed 's/-/\//g')/ $(CUTOFF) result/$*







################################################################
# % = $(Q)/$(chr)/$(ld)
# observed
# jobs/step3/%-fgwas-jobs:
# 	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
# 	@printf "" > $@
# 	@[ -f result/fgwas/$(K)/$*.zscore.gz ] || echo ./make.run-fgwas.R $(shell echo $* | awk -F'/' '{ print ("$(TEMPDIR)/" $$2 "/" $$3) OFS "1KG_EUR/chr" $$2 OFS $(K) OFS "FALSE" OFS ("result/conf/" $$1 "/$(CK)/loco_" $$2 ".txt.gz") OFS ("phenotypes/ukbb_pheno_" $$1 ".txt") OFS ("result/fgwas/$(K)/$*" ) }') >> $@
# 	@[ -f result/fgwas/$(K)/$*.zscore.gz ] || echo ./make.run-fgwas.R $(shell echo $* | awk -F'/' '{ print ("$(TEMPDIR)/" $$2 "/" $$3) OFS "1KG_EUR/chr" $$2 OFS $(K) OFS "TRUE" OFS ("result/conf/" $$1 "/$(CK)/loco_" $$2 ".txt.gz") OFS ("phenotypes/ukbb_pheno_" $$1 ".txt") OFS ("result/fgwas_nn/$(K)/$*" ) }') >> $@

# % = $(chr)/$(ld)
# null
## jobs/step3/%-null-jobs:
##	@for kk in $(K); do [ -f result/null-fgwas_nn/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas-null.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} TRUE result/null-fgwas_nn/$${kk}/$* >> $@; done
##	@for kk in $(K); do [ -f result/null-fgwas/$${kk}/$*.zscore.gz ] || echo ./make.run-fgwas-null.R $(TEMPDIR)/$* 1KG_EUR/chr$(shell echo $* | awk -F'/' '{ print $$1 }') $${kk} FALSE result/null-fgwas/$${kk}/$* >> $@; done

