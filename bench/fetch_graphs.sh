#!/bin/bash

echo "fetch set of common graphs"
cd data
wget -r -nd --no-parent -A '*-mis.tar.gz' "http://sites.nlsde.buaa.edu.cn/~kexu/benchmarks/"

echo "extract graphs"
for file in $(find . -name "*.tar.gz"); do
  tar -xf "$file"
  rm "$file"
done
