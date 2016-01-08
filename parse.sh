#!/bin/bash
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./parse-script.R > tmp_parse.txt 2>&1 &

