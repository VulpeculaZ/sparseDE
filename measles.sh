#!/bin/bash
for ((i = 0; i < 15; i++))
do
    /usr/bin/nohup /bin/nice -n 10 /group/statsoft/R-patched/build-MKL-seq/bin/Rscript --vanilla ./measles-scr04.R $i > tmp_measles_$i.txt 2>&1 &
done
