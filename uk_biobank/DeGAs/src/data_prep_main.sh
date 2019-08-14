#!/bin/bash
set -beEu -o pipefail

if [ $# -lt 3 ] ; then echo "usage: $0 <input dir (root) [ example: /home/scidb/R_code/results/ ] > <output dir (root) [ example: /opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser/static/decomposition ]> <dataset> [loading]" >&2 ; exit 1 ; fi

in_dir=$1
out_dir=$2
dataset=$3
if [ $# -gt 3 ] && [ $4 == "loading" ] ; then loading="loading-true" ; else loading="loading-false" ; fi

variant_tsv=/home/ytanigaw/repos/rivas-lab/decomposition/private_data/variant_and_gene_labels.tsv.gz
threshold=0.01

json_py=$(dirname $0)/data_prep_for_gbe.py

if [ $loading == "loading-true" ] ; then
    python2 $json_py -i $in_dir -o $out_dir -n $dataset -a $variant_tsv -t $threshold -p -j -l
else 
    python2 $json_py -i $in_dir -o $out_dir -n $dataset -a $variant_tsv -t $threshold -p -j
fi

