#!/bin/sh
#
# $Id: response.sh 327 2012-07-13 11:04:42Z jbao $
#
#$ -cwd
#$ -S /bin/bash
# -m eas
# -M jie.bao@frias.uni-freiburg.de
#$ -pe Mmpi 1
#$ -j y
#$ -l h_vmem=500M
#$ -l h_cpu=23:00:00
#$ -R y

#cd $HOME/tool/R-2.15.1/bin
net=scalefree2
size=1000
wd=data/DREAM/gnw/$net/gnw/Size1000/norm_perturbation/mds/
#wd=data/pc12/sohse/120712/
if [ ! -d /export/work/jbao/$wd ]; then
    mkdir -p /export/work/jbao/$wd
fi
for i in `seq 1 100`; do
    for j in `seq 1 50`; do
        filename=$HOME/$wd/pca_${net}-${i}_perturbation-${j}_${size}_normed_euc 
        #filename=$HOME/$wd/pca_${net}_${size}_normed_euc
        prefix=pval_${net}-${i}_perturbation-${j}_$size
        #prefix=pval_${net}_${size}
        # $HOME/tool/R-2.15.1/bin/R CMD BATCH --no-timing "--args filename='$filename' outdir='/export/work/jbao/$wd' prefix='$prefix'" $HOME/github/mds/R/response.R $HOME/$wd/Rout.$i.$SGE_TASK_ID
        R CMD BATCH --no-timing "--args filename='$filename' outdir='/export/work/jbao/$wd' prefix='$prefix'" $HOME/github/mds/R/response.R $HOME/$wd/Rout.$i.$j
        mv /export/work/jbao/$wd/${prefix}.dat $HOME/$wd
    done
done
