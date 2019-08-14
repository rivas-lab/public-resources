#!/bin/bash

phenotype_dir=$1
ontology=$2
score="BFold"

while read d ; do Rscript circular_bar_great.R $d $ontology $score ; done < <( find $phenotype_dir -type d -maxdepth 1 | egrep -v "${phenotype_dir}\$")


