#!/bin/sh
#
# $Id: mds.sh 327 2012-07-13 11:04:42Z jbao $
#
#$ -cwd
#$ -S /bin/bash
#$ -pe Mmpi 8
#$ -l h_vmem=500M
#$ -l h_cpu=23:00:00
#$ -R y
#$ -j y

net=$1
size=1000
wd="$HOME/data/DREAM/gnw/$net/gnw/Size1000/test/mds"
#wd="$HOME/data/pc12/sohse/120712/"
#wd="$HOME/work/DREAM/DREAM3_in_silico_challenge/Size$size/gnw/mds"
#fname="${net}_${size}_rewiring_normed.dat"
#core=`echo $fname | awk -F. '{ print $1; }'`

#interval=200
#arg=$1

cd $HOME/github/mds
i=$SGE_TASK_ID
#i=1

#: <<'END'
for idx in `seq 1 1`; do
#for i in `seq $((($1-1)*$interval+1)) $(($1*$interval))`;do
fname=`ls ${wd}/${net}-${idx}_perturbation-${i}_${size}_normed.dat`
#fname=`ls ${wd}/${net}-full_perturbation-${i}_${size}_normed.dat`
#fname="${wd}/${net}_normed_window$i.dat"
#echo $fname
core=`echo $fname | awk -F. '{ print $1; }'`
core=`echo $core | awk -F/ '{ print $12; }'`
#size=`echo $fname | awk -F_ '{ print $2; }'`
# $HOME/tool/openmpi-1.3.2/bin/mpirun -np $NSLOTS 
$HOME/tool/openmpi-1.3.2/bin/mpirun -n $NSLOTS ./dtw $fname 1 /export/work/jbao/dist_matrix.bin.$net.$i
# $HOME/tool/openmpi-1.3.2/bin/mpirun --mca pls_gridengine_verbose 1 -np $NSLOTS ./dtw $fname 1 /export/work/jbao/dist_matrix.bin.$net.$i
if [ $? -ne 0 ]; then
    echo 'File not found, now exiting...'
    exit 1
fi
mv /export/work/jbao/dist_matrix.bin.$net.$i $HOME/tmp/
#echo $core
time ./hitmds2 100 2 0.04 0 $HOME/tmp/dist_matrix.bin.$net.$i $size $wd/mds_${core}_euc
./pca $wd/mds_${core}_euc $size 2 $wd/pca_${core}_euc
rm $HOME/tmp/dist_matrix.bin.$net.$i
done
#mail -s "Job done - MDS random 1000" jie.bao@gmail.com < /dev/null
#END

#i=6
# $HOME/tool/openmpi-1.3.2/bin/mpirun -np $NSLOTS ./dtw $wd/${net}.csv 1 /export/work/jbao/${net}_${size}.bin
#mv /export/work/jbao/${net}_${size}.bin $HOME/tmp/
#./hitmds2 100 2 0.04 0 $HOME/tmp/${net}_${size}.bin $size $wd/mds_${net}_${size}_normed_euc
#./pca $wd/mds_${net}_${size}_normed_euc $size 2 $wd/pca_${net}_${size}_normed_euc
#rm $HOME/tmp/${net}_${size}.bin
