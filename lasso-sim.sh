#!/bin/bash
for ((i = 0; i <= 19; i++))
do
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./simulation/delay-tv-lasso.R $i > tmp_$i.txt 2>&1 &
done
