#!/bin/bash
set -beEuo pipefail

bed_tar=$1

echo "scp"

hoxa_dir=repos/decomposition-great

scp $bed_tar hoxa:$hoxa_dir
ssh hoxa -t ssh atp -t "bash $hoxa_dir/decomposition-great.sh "
scp hoxa:$hoxa_dir/results.tar.gz $(dirname $bed_tar)
cd $(dirname $bed_tar)
tar xzvf results.tar.gz
cd -

