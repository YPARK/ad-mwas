
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

## combine confounder correction and MWAS



