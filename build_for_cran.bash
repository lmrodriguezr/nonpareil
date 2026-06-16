#!/bin/bash

rm -f Nonpareil_[34].*.tar.gz
cd $(dirname $0)
( cd utils && ./build.bash )
R CMD build utils/Nonpareil
R CMD check --as-cran Nonpareil_[34].*.tar.gz

