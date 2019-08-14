#!/bin/bash

dir=/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser/static/decomposition

tmpfile=$(mktemp)

find $dir -maxdepth 1 -type d | sed -s "s%$dir%%g"| cut -c2- | sort | awk 'length($0) > 0' > $tmpfile
sudo mv $tmpfile $dir/decomposition_datasets.lst
