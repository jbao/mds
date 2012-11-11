// $Id: dtw.cpp 192 2010-10-06 22:11:34Z jbao $

#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
//#include <gsl/gsl_matrix.h>
#include "dtw.h"

DTW::DTW(vector<double>& s, vector<double>& t) : s_(s), t_(t){
    //return_ = new double();
}

DTW::~DTW(){
    //delete return_;
}

double DTW::distance(){
    double gap = 0;

    // actual cost at each position
    const int row = s_.size();
    const int col = t_.size();
    double cost[row][col];
    int i = 0, j;
    for (; i < s_.size(); ++i)
        for (j = 0; j < t_.size(); ++j)
            cost[i][j] = fabs(s_[i]-t_[j]);

    // distance matrix for backtracing
    double dist[row][col];

    // initialize the 1st row and column of dist
    dist[0][0] = cost[0][0];
    double dist_max = dist[0][0];
    for (i = 0; i < s_.size(); ++i){
        dist[i][0] = cost[i][0] + dist[i-1][0] + gap;
        if (dist[i][0] > dist_max)
            dist_max = dist[i][0];
    }
    for (j = 0; j < t_.size(); ++j){
        dist[0][j] = cost[0][j] + dist[0][j-1] + gap;
        if (dist[0][j] > dist_max)
            dist_max = dist[0][j];
    }

    // fill the distance matrix
    double value;
    for (i = 1; i < s_.size(); ++i){
        for (j = 1; j < t_.size(); ++j){
            value = ((dist[i-1][j] + gap) < (dist[i][j-1] + gap)) ? dist[i-1][j] :
                dist[i][j-1];
            value = (value < dist[i-1][j-1]) ? value : dist[i-1][j-1];
            value += cost[i][j];
            dist[i][j] = value;
            if (dist[i][j] > dist_max)
                dist_max = dist[i][j];
            //delete value;
        }
    }
    
    double return_ = dist[row - 1][col - 1]; /// dist_max; 
    /*
    std::cerr << row << ' ' << col << std::endl;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) 
            std::cerr << cost[i][j] << ' ';
        std::cerr << std::endl;
    }
    */
    //delete cost, dist;

    return return_;

}
