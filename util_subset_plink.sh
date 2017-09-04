#!/bin/bash -l

plink_hdr=$1
chr=$2
start_pos=$3
end_pos=$4
out_hdr=$5

mkdir -p $(dirname $out_hdr)

./bin/plink --bfile $plink_hdr --make-bed --chr $chr \
	    --from-bp $start_pos --to-bp $end_pos \
	    --out $out_hdr
