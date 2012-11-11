#!/bin/sh
#
# $Id: responseBatch.sh 271 2011-10-20 14:06:08Z jbao $

# determine the number of jobs
filename='/home/jbao/data/hdf/pca_hdf_nhk_25456_normed_euc'
wd=data/hdf/
line=`awk 'END { print NR }' $filename`
qsub -t 1 response.sh $filename $wd
#for i in `seq 1 $line`;do
#    if [ ! -f $HOME/$wd/pval_$i.Rdat ];then
#        qsub -t $i response.sh $filename $wd
#    fi
#done

#filename='/home/jbao/data/data_aceton_w_mds.txt'
#wd=data/rage/
#line=`awk 'END { print NR }' $filename`
#for i in `seq 1 $line`;do
#    if [ ! -f $HOME/$wd/pval_$i.Rdat ];then
#        qsub -t $i response.sh $filename $wd
#    fi
#done

#filename='/home/jbao/data/cfu/mds_cfu_epo_21464_euc'
#wd=data/cfu/
#line=`awk 'END { print NR }' $filename`
#for i in `seq 1 $line`;do
#    if [ ! -f $HOME/$wd/pval_$i.Rdat ];then
#        qsub -t $i response.sh $filename $wd
#    fi
#done

#filename='/home/jbao/data/yeast/pca_yeast_alpha_6178_normed_euc'
#wd=data/yeast/
#line=`awk 'END { print NR }' $filename`
#for i in `seq 1 $line`;do
#    if [ ! -f $HOME/$wd/pval_$i.Rdat ];then
#        qsub -t $i response.sh $filename $wd
#    fi
#done
