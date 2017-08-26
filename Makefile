
## De Jager et al. 2014
DDIR := /broad/dejagerlab/cogdec/datasets/ROSMAP/DNA_Methylation/DNA_Methylation_Brain_DLPFC_450K/frozenDataSet

RAW := /broad/compbio/eaton/dnam_data/full_data/FinalReport_all_fields_768.txt.gz
IDs := /broad/compbio/eaton/dnam_data/full_data/master_ids.FromLori20121102.csv.xz
PHENO := phenotype/pheno_cov_n3033_032315.csv
GENO := genotype/step08/BED

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
 $(foreach chr, $(CHR), data/raw/chr$(chr)-probes.txt.gz) \
 $(foreach chr, $(CHR), data/raw/chr$(chr)-data.txt.gz) \
 $(foreach chr, $(CHR), data/meth/chr$(chr)-logit.ft) \
 $(foreach chr, $(CHR), data/geno/chr$(chr).geno.mat.ft) \
 data/raw/samples.txt data/raw/masterid.csv

data/raw/all-probes.txt.gz:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $(RAW) | awk -F'\t' '{ print $$2 FS NR }' | gzip > $@

data/raw/chr%-probes.txt.gz: data/raw/all-probes.txt.gz data/probes/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make.chr.probes.R $^ $@

data/raw/samples.txt:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $(RAW) | head -n1 | tr '\t' '\n' | awk '/AVG_Beta/ { gsub(/.AVG_Beta/,""); print }' > $@

data/raw/masterid.csv:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(IDs) | xz -d > $@

data/raw/matched.samples.txt: data/raw/samples.txt data/raw/masterid.csv
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make.matched.samples.R $^ $@

data/raw/chr%-data.txt.gz: data/raw/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	zcat $< | awk 'NR > 1 { printf "," } { printf $$NF }' > data/raw/temp-chr$*.rows
	zcat $(RAW) | awk -F'\t' -vROWSF=data/raw/temp-chr$*.rows -f util_subset_final_report.awk | awk -F'\t' -f util_transpose.awk | gzip -c > $@
	rm data/raw/temp-chr$*.rows

## feather file for faster access
data/meth/chr%-logit.ft: data/meth/chr%-data.ft
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./util.take.logit.R $< $@

data/meth/chr%-data.ft: data/raw/chr%-data.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./util.txt2feather.R $< $@

data/geno/chr%.geno.mat.ft: data/raw/matched.samples.txt
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./util.plink2feather.R genotype/step08/BED/chr$* $< data/geno/chr$*


################################################################
## Generate QTL data in temporary directory
step3: jobs/qtl-data-jobs.txt.gz

clear-step3:
	[ -f jobs/qtl-data-jobs.txt.gz ] && rm jobs/qtl-data-jobs.txt.gz

TEMPDIR := /broad/hptmp/ypp/AD/mwas/

$(TEMPDIR):
	[ -d $@ ] || mkdir -p $@

## Create methylation job lists
jobs/qtl-data-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-data-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N qtl.data -binding "linear:1" -q short -l h_vmem=8g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

CHUNK := 20 # group every ~20 CpGs as one job
CTRL := 5   # control CpGs for each CpG
jobs/temp-qtl-data-%-jobs.txt.gz: data/probes/chr%-probes.txt.gz $(TEMPDIR)
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make_job_segments.awk -vNTOT=$$(zcat $< | wc -l) -vCHUNK=$(CHUNK) | awk '{ print "./data.qtl.R" FS $* FS $$1 FS $(CTRL) FS ("$(TEMPDIR)/$*/" NR "-data") }' | gzip > $@


################################################################
## Confounder correction and QTL calling
step4: jobs/qtl-run-jobs.txt.gz jobs/qtl-perm-jobs.txt.gz

clear-step4:
	[ -f jobs/qtl-run-jobs.txt.gz ] && rm jobs/qtl-run-jobs.txt.gz
	[ -f jobs/qtl-perm-jobs.txt.gz ] && rm jobs/qtl-perm-jobs.txt.gz

jobs/qtl-run-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-run-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N qtl.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

jobs/temp-qtl-run-%-jobs.txt.gz: data/probes/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make_job_segments.awk -vNTOT=$$(zcat $< | wc -l) -vCHUNK=$(CHUNK) | awk '{ print "./run.qtl.R" FS ("$(TEMPDIR)/$*/" NR "-data") FS ("result/qtl/chr$*/b" NR "/qtl") }' | gzip > $@

jobs/qtl-perm-jobs.txt.gz: $(foreach chr, $(CHR), jobs/temp-qtl-perm-$(chr)-jobs.txt.gz)
	cat $^ > $@
	@[ $$(zcat $@ | wc -l) -lt 1 ] || qsub -t 1-$$(zcat $@ | wc -l) -N qtl.data -binding "linear:1" -q short -l h_vmem=4g -P compbio_lab -V -cwd -o /dev/null -b y -j y ./run_rscript.sh $@
	rm $^

jobs/temp-qtl-perm-%-jobs.txt.gz: data/probes/chr%-probes.txt.gz
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	./make_job_segments.awk -vNTOT=$$(zcat $< | wc -l) -vCHUNK=$(CHUNK) | awk '{ print "./run.qtl.R" FS ("$(TEMPDIR)/$*/" NR "-data") FS ("result/perm/chr$*/b" NR "/perm") FS "TRUE" }' | gzip > $@
