#!/bin/bash
for ((i = 1; i <= 50; i++))
do
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./simulation/sim.tv.fused02.R $i > tmp_tv$i.txt 2>&1 &
done
