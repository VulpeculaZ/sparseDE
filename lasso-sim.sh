#!/bin/bash
for ((i = 1; i <= 50; i++))
do
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./simulation/sim.2d.lasso.R $i > tmp_$i.txt 2>&1
done