// $Id: dtw.h 118 2010-01-27 10:28:14Z jbao $

#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <fstream>

using std::vector;
using std::ifstream;

typedef vector< vector<double> > DATA;

void load_dat(ifstream& file, DATA& data, vector<double>& time_points);
void load_mat(ifstream& file, DATA& data);
double getEuclideanDistance(vector<double>& s, vector<double>& t);
double getEuclideanDistance(double s, double t);

#endif
