mds
===

Multidimensional scaling

Adapted from HiT-MDS-2 http://pgrc.ipk-gatersleben.de/seeds/analysis_tools.php

Details
=======

A customized program (`dtw`) is first used to calculate the Euclidean distances
between time series, the code is parallelized by MPI, but one only gets a 
significant speed-up when dealing with large networks.

`hitmds2` is then applied on the output distance matrix and generates a 
two-dimensional projection of individual entities.

`pca` then rotates the projection such that the first dimension represents the
largest variation in data.

To run in parallel,
        
    for i in `seq 1 100`;do qsub -t 1-50 mds.sh scalefree2 $i;done

The array job index 1-50 denote different perturbations and the loop index 1-100 
denote different topologies.
