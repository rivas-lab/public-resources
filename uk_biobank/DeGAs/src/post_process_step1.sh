#!/bin/bash
set -beEuo pipefail

npz_file=$1
#npz_file=../private_data/npz/dev_PTVsNonMHC_z_nonCenter_p0001_100PCs.npz

python post_process_step1.py -i $npz_file -o ../private_data/results/ -p ../public_data/phenotype_of_interest.lst
