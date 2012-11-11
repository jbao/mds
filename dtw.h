// $Id: dtw.h 192 2010-10-06 22:11:34Z jbao $

#ifndef DTW_H
#define DTW_H

#include <vector>
//#include <gsl/gsl_matrix.h>

using std::vector;

class DTW{

public:
    DTW(vector<double>& s, vector<double>& t);
    ~DTW();
    double distance();

private:
    vector<double> s_, t_;
    //double return_;

};

#endif
