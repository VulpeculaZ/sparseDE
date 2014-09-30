#!/bin/bash
for ((i = 0; i < 24; i++))
do
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./simulation/measles-scr.R $i > tmp_measles_$i.txt 2>&1 &
done
