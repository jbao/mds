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
#$ -l hostname=head

#cd $HOME/tool/R-2.15.1/bin
export LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/lib/jvm/java-6-sun/lib/amd64/server:/usr/lib/jvm/java-6-sun/lib/amd64:/usr/lib/x86_64-linux-gnu/jni:/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu:/home/jbao/tool/graphviz-2.26.3/lib:/home/jbao/tool/openmpi-1.3.2/lib:/home/jbao/tool/rpy-2.1.4/lib:/home/jbao/tool/libsbml-4.2.0/lib:/home/jbao/tool/getopt_pp:/home/jbao/tool/armadillo-2.4.3/:/home/sweber/.local/lib:/usr/lib
echo $LD_LIBRARY_PATH
net=$1
size=1000
wd=data/DREAM/gnw/$net/gnw/Size1000/norm_perturbation/mds/
#wd=data/pc12/sohse/120712/
if [ ! -d /export/work/jbao/$wd ]; then
    mkdir -p /export/work/jbao/$wd
fi
#for i in `seq 1 1`; do
i=$2
filename=$HOME/$wd/pca_${net}-${i}_perturbation-${SGE_TASK_ID}_${size}_normed_euc 
#filename=$HOME/$wd/pca_${net}_${size}_normed_euc
prefix=pval_${net}-${i}_perturbation-${SGE_TASK_ID}_$size
#prefix=pval_${net}_${size}
# $HOME/tool/R-2.15.1/bin/R CMD BATCH --no-timing "--args filename='$filename' outdir='/export/work/jbao/$wd' prefix='$prefix'" $HOME/github/mds/R/response.R $HOME/$wd/Rout.$i.$SGE_TASK_ID
R CMD BATCH --no-timing --verbose "--args filename='$filename' outdir='/export/work/jbao/$wd' prefix='$prefix'" $HOME/github/mds/R/response.R $HOME/$wd/Rout.$i.$SGE_TASK_ID
mv /export/work/jbao/$wd/${prefix}.dat $HOME/$wd
#done
