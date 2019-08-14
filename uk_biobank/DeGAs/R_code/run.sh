#!/bin/bash

# you can loop over different conditions (bim_filters, etc.)

#skipR="skipR"
skipR="runR"
algorithm="tsvd"
R_script_base="tsvd.R"

for z_or_lor in z ; do
for center_tsvd in center ; do
for num_components in 90 100 110 ; do
for pvalue in 0.001 ; do
for bim_filter in PTVsNonMHC allNonMHC ; do
  bash $(dirname $(readlink -f $0))/run-sub.sh ${bim_filter} ${z_or_lor} ${center_tsvd} ${pvalue} ${num_components} ${skipR} ${algorithm} ${R_script_base}
done
done
done
done
done

