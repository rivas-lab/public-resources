#!/bin/bash
set -beEu -o pipefail

if [ $# -lt 1 ] ; then echo "usage: $0 [bim_filter: codingNonMHC] [z (default) | lor] [nonCenter (default) | non-center] [p-val: 0.001] [num components: 100] [skipR] [method: tsvd (default) | ssvd] [script name: tsvd.R] [value_filter_mode: default, pvalue, or se] [additional suffix]" >&2 ; exit 1 ; fi

if [ $# -gt 0 ] ; then bim_filter=$1     ; else bim_filter="codingNonMHC" ; fi
if [ $# -gt 1 ] ; then z_or_lor=$2       ; else z_or_lor="z" ; fi
if [ $# -gt 2 ] ; then center_tsvd=$3    ; else center_tsvd="nonCenter" ; fi
if [ $# -gt 3 ] ; then pvalue=$4         ; else pvalue="0.001" ; fi
if [ $# -gt 4 ] ; then num_components=$5 ; else num_components=100 ; fi
if [ $# -gt 5 ] && [ $6 == "skipR" ] ; then skipR=true ; else skipR=false ; fi
if [ $# -gt 6 ] && [ $7 == "ssvd" ] ; then algorithm=ssvd ; else algorithm=tsvd ; fi
if [ $# -gt 7 ] ; then Rscript_base=$8    ; else Rscript_base="tsvd.R" ; fi
if [ $# -gt 8 ] ; then value_filter_mode=$9 ; else value_filter_mode="default" ; fi
if [ $# -gt 9 ] ; then custom_suffix=$10 ; else custom_suffix="" ; fi

Rscript=$(readlink -e "$(dirname $(readlink -e $0))/${Rscript_base}")

dataset="dev_${bim_filter}_${z_or_lor}_${center_tsvd}_p$(echo $pvalue | sed -e 's/\.//g')_${algorithm}_${num_components}PCs${custom_suffix}_$(date +%Y%m%d)"

echo "dataset: ${dataset}"

out_d=/home/scidb/R_code/results/${dataset}
gbe_d=/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser/static/decomposition

# create output directory (with logs subdir)
if [ ! -d $out_d/logs ] ; then mkdir -p $out_d/logs ; fi
timestamp_f=$out_d/metadata.txt
log_file=$out_d/logs/tsvd.log
echo "start time: $(date)" | tee $timestamp_f
cp $0       $out_d/logs/
cp $Rscript $out_d/logs/

# dump git log
echo "$ git log -n1"
cd $(dirname $(readlink -e $0))
echo "script git version: $(git log -n1| awk 'length($0)>0'| tr "\n" "\t")" | tee -a $timestamp_f
cd -

if [ $skipR != "true" ] ; then
  # TSVD computation
  echo "running $ Rscript $Rscript $dataset $bim_filter $z_or_lor $center_tsvd $pvalue $num_components $algorithm $value_filter_mode ... "
  cd $(dirname $Rscript)
  Rscript $Rscript $dataset $bim_filter $z_or_lor $center_tsvd $pvalue $num_components $algorithm $value_filter_mode 2>&1 | tee -a $log_file
  cd -
  echo "TSVD computation finished on $(date)"  | tee -a $timestamp_f
else
  echo "skipR flag is enabled. skipping TSVD computation and applying post-processing to the existing results file from TSVD ..."
fi

if [ ! -e $out_d/ap_variant_idx.tsv ] || [ ! -e $out_d/ap_icd_idx.tsv ] || [ ! -e $out_d/total_inertia.txt ] ; then
  # Export variant and phenotype labels
  iquery -o tsv+:l -aq "scan(ap_variant_idx_${dataset})"  | gzip -9 > $out_d/ap_variant_idx.tsv.gz
  iquery -o tsv+:l -aq "scan(ap_icd_idx_${dataset})"      | gzip -9 > $out_d/ap_icd_idx.tsv.gz
  # Compute and export total inertia 
  iquery -o tsv    -aq "aggregate(apply(ap_icd_var_matrix_${dataset}, z_sq, z*z), sum(z_sq))" > $out_d/total_inertia.txt

  # Export correlation (loading) 
  #  - ap_icd_svd_cor1 - [ICD by PC]
  #  - ap_icd_svd_cor2 - [VARIANT by PC]
  # iquery -o tsv+:l -aq "scan(ap_icd_cor1_${dataset})" | gzip -9 > $out_d/ap_icd_svd_cor1.tsv.gz
  # iquery -o tsv+:l -aq "scan(ap_icd_cor2_${dataset})" | gzip -9 > $out_d/ap_icd_svd_cor2.tsv.gz
fi

if [ $algorithm == "tsvd" ] && [ ! -e $out_d/ap_icd_var_tsvd.tsv ]  ; then
  # Export TSVD decomposed matrix
  iquery -o tsv+:l -aq "scan(ap_icd_var_tsvd_${dataset})" | gzip -9 > $out_d/ap_icd_var_tsvd.tsv.gz
fi

# export to python numpy file and JSON files
sudo bash $(dirname $(dirname $(readlink -e $0)))/src/data_prep_main.sh $(dirname $out_d) $gbe_d $dataset # loading

cd $gbe_d
sudo tar czf ${dataset}.tar.gz ${dataset}
cd -

echo $gbe_d/${dataset}.tar.gz
echo ""

