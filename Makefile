
## De Jager et al. 2014
DDIR := /broad/dejagerlab/cogdec/datasets/ROSMAP/DNA_Methylation/DNA_Methylation_Brain_DLPFC_450K/frozenDataSet

RAW := /broad/compbio/eaton/dnam_data/full_data/FinalReport_all_fields_768.txt.gz
IDs := /broad/compbio/eaton/dnam_data/full_data/master_ids.FromLori20121102.csv.xz
PHENO := phenotype/pheno_cov_n3033_032315.csv
GENO := genotype/step08/BED
CPROBE := /broad/compbio/eaton/dnam_data/full_data/DLPFC_450KMethy_FinalDataSet_8.2.11/Control\ Probe\ Profile.txt.xz

PLINKZIP := https://www.cog-genomics.org/static/bin/plink170815/plink_linux_x86_64.zip

CHR := $(shell seq 1 22)

################################################################
## Copy list of CpGs Q/C'ed by De Jager et al. 2014

step1: $(foreach chr, $(CHR), data/probes/chr$(chr)-probes.txt.gz)

data/probes/chr%-probes.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(DDIR)/ill450kAnno_finalQC_chr$*.txt | tail -n+2 |awk -F'\t' '{ print $$1 FS $$3 FS $$4 }' | sort -k3n | gzip -c > $@

################################################################
## Break down large Matt Eaton's final report file into small pieces
## and save only samples with geneotypes

step2: data/raw/all-probes.txt.gz \
 data/raw/control-probes.txt.gz \
 data/raw/samples.txt data/raw/masterid.csv \
 data/raw/matched.samples.txt.gz \
 data/meth/control.ft \
 $(foreach chr, $(CHR), data/raw/chr$(chr)-probes.txt.gz) \
 $(foreach chr, $(CHR), data/raw/chr$(chr)-beta.txt.gz) \
 $(foreach chr, $(CHR), data/meth/chr$(chr)-logit.ft) \
 data/meth/PC.txt.gz

data/raw/all-probes.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $(RAW) | awk -F'\t' '{ print $$2 FS NR }' | gzip > $@

data/raw/chr%-probes.txt.gz: data/raw/all-probes.txt.gz data/probes/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make.chr.probes.R $^ $@

data/raw/samples.txt:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $(RAW) | head -n1 | tr '\t' '\n' | awk -F'\t' '/AVG_Beta/ { gsub(/.AVG_Beta/,""); print $$0 FS (++id) FS NR }' > $@

data/raw/masterid.csv:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(IDs) | xz -d > $@

data/raw/matched.samples.txt.gz: data/raw/samples.txt data/raw/masterid.csv
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make.matched.samples.R $^ $@

data/raw/chr%-beta.txt.gz: data/raw/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk 'NR > 1 { printf "," } { printf $$NF }' > data/raw/temp-chr$*.rows
	zcat $(RAW) | awk -F'\t' -vROWSF=data/raw/temp-chr$*.rows -f util_subset_final_report.awk | awk -F'\t' -f util_transpose.awk | gzip -c > $@
	rm data/raw/temp-chr$*.rows

## control probes
data/raw/control-probes.txt.gz:
	cat $(CPROBE) | xz -d | gzip > $@

data/meth/control.ft: data/raw/matched.samples.txt.gz data/raw/control-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make.matched.control.R $^ $@

## PCs of most variable CpG probes
data/meth/PC.txt.gz: data/raw/matched.samples.txt.gz $(foreach chr, $(CHR), data/raw/chr$(chr)-beta.txt.gz) $(foreach chr, $(CHR), data/meth/chr$(chr)-logit.ft)
	./make.meth.pc.R

## feather file for faster access
data/meth/chr%-logit.ft: data/raw/chr%-probes.txt.gz data/raw/chr%-beta.txt.gz data/raw/samples.txt
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./util.take.logit.R $^ $@

################################################################
## Generate QTL data in temporary directory
step3: bin/plink \
  $(foreach chr, $(CHR), jobs/segments/$(chr)-jobs.txt) \
  jobs/qtl-data-jobs.txt.gz

clear-step3:
	[ -f jobs/qtl-data-jobs.txt.gz ] && rm jobs/qtl-data-jobs.txt.gz

TEMPDIR := /broad/hptmp/ypp/AD/mwas/qtl/

$(TEMPDIR):
	[ -d $@ ] || mkdir -p $@

jobs/segments/%-jobs.txt: data/probes/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./run.sh ./make.job.segments.R $< > $@

## Create methylation job lists
jobs/qtl-data-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-data-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N qtl.data -binding "linear:1" -q short -l h_vmem=6g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

NCTRL := 3   # 3 control CpGs for each CpG
DATA_QTL := ./make.data.qtl.R
GENO_HDR := ./rosmap-geno/gen/impute/rosmap1709-chr

## % = $(chr)
jobs/temp-qtl-data-%-jobs.txt.gz: jobs/segments/%-jobs.txt data/meth/chr%-logit.ft $(TEMPDIR)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk '{ print "$(DATA_QTL)" FS "data/meth/chr$*-logit.ft" FS "data/raw/chr$*-probes.txt.gz" FS $$1 FS $(NCTRL) FS "$(GENO_HDR)$*" FS ("$(TEMPDIR)/$*/" NR "-data") }' | gzip > $@

################################################################
## Confounder correction and QTL calling
step4: jobs/qtl-run-jobs.txt.gz

clear-step4:
	[ -f jobs/qtl-run-jobs.txt.gz ] && rm jobs/qtl-run-jobs.txt.gz

jobs/qtl-run-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-run-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N qtl.run -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

jobs/temp-qtl-run-%-jobs.txt.gz: jobs/segments/%-jobs.txt
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk '{ print "./make.summary.qtl.R" FS ("$(TEMPDIR)/$*/" NR "-data") FS ("result/qtl/$*/b" NR) FS $$1 }' | gzip > $@

################################################################
## Hi-C guided QTL calling
step5: jobs/hic-qtl-data-jobs.txt.gz

TEMPDIR-HIC := /broad/hptmp/ypp/AD/mwas/hic-qtl/

jobs/hic-qtl-data-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-hic-qtl-data-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N hic-qtl.data -binding "linear:1" -l h_rt=1800 -l h_vmem=6g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

## % = $(chr)
jobs/temp-hic-qtl-data-%-jobs.txt.gz: jobs/segments/%-jobs.txt data/meth/chr%-logit.ft $(TEMPDIR-HIC)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $< | awk '{ print "$(DATA_QTL)" FS "data/meth/chr$*-logit.ft" FS "data/raw/chr$*-probes.txt.gz" FS $$1 FS $(NCTRL) FS "$(GENO_HDR)$*" FS 3e6 FS ("$(TEMPDIR-HIC)/$*/" NR "-data") }' | gzip > $@

$(TEMPDIR-HIC):
	[ -d $@ ] || mkdir -p $@

################################################################
step6: jobs/hic-qtl-run-jobs.txt.gz

jobs/hic-qtl-run-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-hic-qtl-run-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N hic-qtl.run -binding "linear:1" -l h_rt=1:00:00 -l h_vmem=6g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

jobs/temp-hic-qtl-run-%-jobs.txt.gz: jobs/segments/%-jobs.txt
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	printf "" | gzip >> $@
	cat $< | awk '{ print "./make.summary.hic-qtl.R" FS ("$(TEMPDIR-HIC)/$*/data-" NR) FS "hic/CO/chr$*.pairs.gz" FS ("result/hic-qtl/CO/$*/b" NR) }' | gzip >> $@
	cat $< | awk '{ print "./make.summary.hic-qtl.R" FS ("$(TEMPDIR-HIC)/$*/data-" NR) FS "hic/HC/chr$*.pairs.gz" FS ("result/hic-qtl/HC/$*/b" NR) }' | gzip >> $@


################################################################
## Utilities
bin/plink:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	curl $(PLINKZIP) -o bin/plink.zip
	unzip bin/plink.zip -d bin/

