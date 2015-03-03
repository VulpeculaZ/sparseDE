#!/bin/bash
for ((i = 0; i < 20; i++))
do
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./simulation/blowfly-sim-batch.R $i > tmp_nnls_$i.txt 2>&1 &
done
