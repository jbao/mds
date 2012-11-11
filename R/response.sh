#!/bin/sh
#
# $Id: response.sh 327 2012-07-13 11:04:42Z jbao $
#
#$ -cwd
#$ -S /bin/bash
# -m eas
# -M jie.bao@frias.uni-freiburg.de
#$ -pe Mmpi 8
#$ -j y
#$ -l h_vmem=500M
#$ -l h_cpu=23:00:00
#$ -R y

net=$1
size=14530
#wd=data/DREAM/gnw/$net/gnw/full/mds/
wd=data/pc12/sohse/120712/
if [ ! -d /export/work/jbao/$wd ]; then
    mkdir -p /export/work/jbao/$wd
fi
for i in `seq 1 1`; do
    #filename=$HOME/$wd/pca_${net}-full_perturbation-${SGE_TASK_ID}_${size}_normed_euc 
    filename=$HOME/$wd/pca_${net}_${size}_normed_euc
    #prefix=pval_${net}-full_perturbation-${SGE_TASK_ID}_$size
    prefix=pval_${net}_${size}
    R CMD BATCH --no-timing "--args filename='$filename' outdir='/export/work/jbao/$wd' prefix='$prefix'" $HOME/MicroarrayAnalysis/trunk/r/response.R $HOME/$wd/Rout.$SGE_TASK_ID
    mv /export/work/jbao/$wd/${prefix}.dat $HOME/$wd
done
