// $Id: main.cpp 327 2012-07-13 11:04:42Z jbao $

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include "mpi.h"
#include "dtw.h"
//#include "mds.h"
#include "main.h"
//#include <gsl/gsl_spline.h>

void load_dat(std::ifstream& data_file, DATA& data, std::vector<double>& time_points){
    if (!data_file.is_open()){
        std::cerr << "Failed to open file" << std::endl;
        return;
    }

    std::string line;
    char delim = '|';

    // time points
    getline(data_file, line);
    std::istringstream ss(line);
    //std::vector<double> time_points;
    std::string entry;
    //getline(ss, entry, '\t');
    while (getline(ss, entry, delim)){
        time_points.push_back(atof(entry.c_str()));    
        //std::cerr << entry << " ";
    }
    
    // data matrix
    while (getline(data_file, line)){
        std::istringstream ls(line);
        std::vector<double> ln;
        //std::string entry;
        getline(ls, entry, delim);
        //std::cerr << entry << " ";
        while (std::getline(ls, entry, delim)){
            ln.push_back(atof(entry.c_str()));
            //std::cerr << entry << " ";
        }
        data.push_back(ln);
        //std::cerr << std::endl;
    }
}

void load_mat(std::ifstream& data_file, DATA& data){
    if (!data_file.is_open()){
        std::cerr << "Failed to open file" << std::endl;
        return;
    }

    std::string line;
    
    // data matrix
    while (getline(data_file, line)){
        std::istringstream ls(line);
        std::vector<double> ln;
        std::string entry;
        //std::cerr << entry << " ";
        while (std::getline(ls, entry, ' ')){
            ln.push_back(atof(entry.c_str()));
            //std::cerr << entry << " ";
        }
        data.push_back(ln);
        //std::cerr << std::endl;
    }
}

double getEuclideanDistance(vector<double>& s, vector<double>& t) {
    double sum = 0;
    vector<double>::iterator iter_s = s.begin();
    vector<double>::iterator iter_t = t.begin();
    for (; iter_s < s.end(); ++iter_s, ++iter_t) {
        sum += (*iter_s - *iter_t) * (*iter_s - *iter_t);
    }
    return sqrt(sum);
}

double getPearsonCorrelation(vector<double>& x, vector<double>& y) {
    const double TINY = 1.0e-20;
    int j;
    double yt, xt, t, df;
    double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;

    assert(x.size() == y.size());
    int n = x.size();
    for (j = 0; j < n; j++) {
        ax += x[j];
        ay += y[j];
    }
    ax /= n;
    ay /= n;
    for (j = 0; j < n; j++) {
        xt = x[j] - ax;
        yt = y[j] - ay;
        sxx += xt * xt;
        syy += yt * yt;
        sxy += xt * yt;
    }
    double r = sxy / (sqrt(sxx * syy) + TINY);
    return r;
}

/**
 * overloading
 */
double getEuclideanDistance(double s, double t) {
    return fabs(s - t);
}

int main(int argc, char *argv[]){

    vector<DATA> all_data;
    DATA all_time_points;
    int row = 0;//, col = atoi(argv[2]);
    //for (int i = 1; i < argc; ++i) {
    DATA data;
    std::ifstream data_file(argv[1]);
    vector<double> time_points;
    //if (std::string(argv[1]).find("dat") != std::string::npos)
        load_dat(data_file, data, time_points);
    //else
    //    load_mat(data_file, data);
    int n_perturb = atoi(argv[2]);
    char *outfile = argv[3];
    row = data.size()/n_perturb;
    //}
    DATA interpol_data;
    double alpha = 0; 
    //int col = data[0].size();    // number of time points
    //std::cerr << row << " " << col << std::endl;

    /*
    // interpolation
    double interpol_step = time_points[1] - time_points[0];
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, col);
    //std::cerr << time_points[0] << std::endl;
    //std::cerr << data[0][0] << std::endl;
    for (int i = 0; i < row; ++i){
        gsl_spline_init(spline, &time_points[0], &data[i][0], col);
        vector<double> gene;
        for (double xi = time_points[0]; xi <= time_points[col-1]; 
            xi += interpol_step){
            
            double yi = gsl_spline_eval(spline, xi, acc);
            //std::cout << xi << " " << yi << std::endl;
            gene.push_back(yi);
        
        }
        interpol_data.push_back(gene);
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    */
    //std::cerr << interpol_data.size() << " " << interpol_data[0].size()
    //    << std::endl;
    //return 0;
    
    //double **data = new double*[row];
    //char filename[50];
    //sprintf(filename, "baf3_treated_%d.dat", row);
    //std::string filename = "baf3_treated_10.dat";
    int i, j, k;
    
    int half = static_cast<int>(row/2);
    //int nScore = nRow[rank]*(row+1);
    //double *all_score = new double[nScore]; 
    //double *buffer = new double[half*(row+1)];

    int rank, size;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > half && rank == 0){
        std::cerr << "Please run with " << half << 
            " processes" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        //return 0;
    }

    // determine the rows to be read by each process
    int **iRow = new int*[size];
    int *nRow = new int[size];
    for (i = 0; i < size; ++i){ // for all processes
        if (i < (row/2)%size)
            nRow[i] = static_cast<int>(row/2/size)+1;
        else
            nRow[i] = static_cast<int>(row/2/size);
        
        iRow[i] = new int[nRow[i]];
        for (j = 0; j < nRow[i]; ++j){ // for all rows
            if (iRow[i][j] == static_cast<int>(row/2))
                break;
            iRow[i][j] = j*size+i%size;
            //std::cerr << iRow[i][j] << " ";
        }
        //std::cerr << std::endl;
    }

    //double **dist;
    //if (rank == 0){
    //    gsl_matrix *dist = gsl_matrix_calloc(row, row);
    //}
    //gsl_matrix_set_zero(dist);
    //vector<double> s, t;
    DTW* dtw_obj;
    double score;

    int recv_counts[size], displs[size];
    for (i = 0; i < size; ++i){
        //recv_counts[i] = nRow[i]*(row+1);
        recv_counts[i] = row * (row / size);
        displs[i] = recv_counts[i] * i;
        //displs[i] = 0;
        //if (i > 0){
        //    for (j = 0; j < i; ++j){
        //        displs[i] += nRow[j]*(row+1);
        //    }
        //}
    }

    MPI_File fh;
    MPI_Offset my_offset = rank * sizeof(double) * (row / size)
        * row;
    MPI_File_open(MPI_COMM_WORLD, outfile, MPI_MODE_CREATE 
        | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_seek(fh, my_offset, MPI_SEEK_SET);

    // DTW
    int iScore = displs[rank];
    for (i = rank*(row/size); i < (rank+1)*(row/size); ++i){ // for all rows 
        //int iScore = displs[rank] + i*(row+1);
        for (j = 0; j < row; ++j){ // for all cols
            //s = vector<double>(col);
            //t = vector<double>(col);
            //std::cerr << iRow[rank][i] << std::endl;
            //for (k = 0; k < col; ++k){ // for all time points
            //    s.push_back(data[iRow[rank][i]][k]);
            //    t.push_back(data[j][k]);
                //std::cerr << data[iRow[rank][i]][k] << " "
                //    << data[j][k] << std::endl;
            //}
            //dtw_obj = new DTW(data[i], data[j]);
            //score = dtw_obj->distance();
            //score = getEuclideanDistance(data[i][col], data[j][col]);
            //score = getEuclideanDistance(data[i], data[j]);
            //score = alpha * (getEuclideanDistance(data[i], data[j]) +
            //    getEuclideanDistance(data[i+row], data[j+row])) + 
            //    (1-alpha) * (getEuclideanDistance(data[i], data[j+row]) +
            //    getEuclideanDistance(data[i+row], data[j])); 
            //score = (getEuclideanDistance(data[i], data[j]) +
            //    getEuclideanDistance(data[i+row], data[j+row]) +
            //    getEuclideanDistance(data[i+2*row], 
            //    data[j+2*row])) / 3;
            score = 0;
            int iPerturb;
            for (iPerturb = 0; iPerturb < n_perturb; ++iPerturb)
                score += getEuclideanDistance(data[i+iPerturb*row],
                        data[j+iPerturb*row]);
                //score += getPearsonCorrelation(data[i+iPerturb*row],
                //        data[j+iPerturb*row]);
            score /= n_perturb;
            //std::cout << "rank " << rank << ": " << i 
            //    << " " << j << " " << std::fixed << score 
            //    << std::endl;
            //buffer[iScore++] = score;
            MPI_File_write(fh, &score, 1, MPI_DOUBLE, &status);
            //delete dtw_obj;
        }
        /*
        // bottom of the distance matrix
        for (j = row-iRow[rank][i]-1; j < row; ++j){ // for all cols
            //s = vector<double>(col);
            //t = vector<double>(col);
            //for (k = 0; k < col; ++k){ // for all time points
            //    s.push_back(data[row-iRow[rank][i]-1][k]);
            //    t.push_back(data[j][k]);
            //}
            dtw_obj = new DTW(data[row-iRow[rank][i]-1], data[j]);
            score = dtw_obj->distance();
            std::cout << "rank " << rank << ": " << row-iRow[rank][i]-1 
                << " " << j << " " << std::fixed << score << std::endl;
            buffer[iScore++] = score;
            //MPI_Send(&score, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
            //MPI_Wait(&request, &status);
            //if (rank == 0){
            //    MPI_Recv(&buffer, 1, MPI_DOUBLE, rank, 123, MPI_COMM_WORLD, 
            //        &status);
            //    dist[j][row-iRow[rank][i]-1] = dist[row-iRow[rank][i]-1][j] 
            //        = buffer;
            //}
            delete dtw_obj;
        }
        */
    }

    // communication
    //MPI_Gatherv(&buffer[displs[rank]], recv_counts[0], 
    //    MPI_DOUBLE, &buffer[0], 
    //    recv_counts, displs, MPI_DOUBLE, 0,
    //    MPI_COMM_WORLD);
    std::cerr << "rank " << rank << ": success" << std::endl;
    if (rank == 0){
        
        //double *dist = new double[row*row];
        //for (i = 0; i < row; ++i){
            //std::cerr << i << std::endl;
            //dist[i] = new double[row];   
        //}

        //gsl_matrix *dist = gsl_matrix_calloc(row, row);
        /*
        // convert buffer to dist
        int iBuffer = 0;
        for (int iRank = 0; iRank < size; ++iRank){
            //std::cerr << sizeof(*buffer)/sizeof(double) << std::endl;
            for (i = 0; i < nRow[iRank]; ++i){ // for all rows 
                int ref = iRow[iRank][i];
                for (j = ref; j < row; ++j){ // for all cols
                    //std::cerr << buffer[iBuffer] << std::endl;
                    //dist[j*row+ref] = dist[ref*row+j] = 
                        //buffer[iBuffer++];
                    gsl_matrix_set(dist, j, ref, buffer[iBuffer++]);
                    gsl_matrix_set(dist, ref, j, buffer[iBuffer++]);
                }
                for (j = row-ref-1; j < row; ++j){ // for all cols
                    //std::cerr << buffer[iBuffer] << std::endl;
                    //dist[j*row+row-ref-1] = 
                        //dist[(row-ref-1)*row+j] = buffer[iBuffer++];
                    gsl_matrix_set(dist, j, row-ref-1,  buffer[iBuffer++]);
                    gsl_matrix_set(dist, row-ref-1, j, buffer[iBuffer++]);
                }
            }
        }
        
        delete[] buffer;

        //}
        */
        // odd number of total genes
        if (row % size != 0){
            
            my_offset = sizeof(double) * (row / size)
                * row * size;
            MPI_File_seek(fh, my_offset, MPI_SEEK_SET);

            //int ref = static_cast<int>(row/2);
            iScore = row * ((row / size) * size);
            for (i = size*(row/size); i < row; ++i){ // for all rows
                for (j = 0; j < row; ++j){ // for all cols
                    //s = vector<double>(col);
                    //t = vector<double>(col);
                    //for (k = 0; k < col; ++k){ // for all time points
                    //    s.push_back(data[ref][k]);
                    //    t.push_back(data[j][k]);
                    //}
                    //dtw_obj = new DTW(data[i], data[j]);
                    //score = dtw_obj->distance();
                    //score = getEuclideanDistance(data[i][col], data[j][col]);
                    //score = getEuclideanDistance(data[i], data[j]);
                    //score = alpha * (getEuclideanDistance(data[i], data[j]) +
                    //    getEuclideanDistance(data[i+row], data[j+row])) + 
                    //    (1-alpha) * (getEuclideanDistance(data[i], data[j+row]) +
                    //    getEuclideanDistance(data[i+row], data[j])); 
                    //score = (getEuclideanDistance(data[i], data[j]) +
                    //    getEuclideanDistance(data[i+row], data[j+row]) +
                    //    getEuclideanDistance(data[i+2*row], 
                    //    data[j+2*row])) / 3;
                    score = 0;
                    int iPerturb;
                    for (iPerturb = 0; iPerturb < n_perturb; ++iPerturb)
                        score += getEuclideanDistance(data[i+iPerturb*row],
                                data[j+iPerturb*row]);
                        //score += getPearsonCorrelation(data[i+iPerturb*row],
                        //        data[j+iPerturb*row]);
                    score /= n_perturb;
                    //std::cout << "rank " << rank << ": " << i 
                    //    << " " << j << " " << std::fixed << score << std::endl;
                    //buffer[iScore] = score;
                    MPI_File_write(fh, &score, 1, MPI_DOUBLE, &status);
                    //delete dtw_obj;
                }
                    
            }
        }
    
    //for (i = 0; i < row; ++i)
    //    for (j = 0; j < i; ++j)
    //        dist[i][j] = dist[j][i];


    //gsl_matrix_free(data);
    //gsl_matrix_free(dist);

    // write to file
    //if (rank == 0){
        //std::ofstream out_file("dist_matrix.bin", std::ios::out | 
	        //std::ios::binary);
	//double *dist_data = new double[row*row];
	//out_file << "#" << -row << " " << -col << std::endl; // for hit-mds
        //for (i = 0; i < row; ++i){
	        //out_file.write((char*)&dist[i], sizeof(dist[i]));
            //for (j = 0; j < row; ++j){
                //out_file << std::fixed << dist[i][j] << " ";
		//dist_data[i*row+j] = dist[i][j];
            //}    
            //out_file << std::endl;
        //}
        //out_file << -2 << std::endl; // target dim, for hit-mds
        //gsl_matrix_view m = gsl_matrix_view_array(buffer, row, row);
        //FILE *f = fopen("dist_matrix.bin", "wb");
        //gsl_matrix_fwrite(f, &m.matrix);
        //fclose(f);
	//out_file.close();
	//gsl_matrix_free(dist);
    }

    MPI_File_close(&fh);

    // free memory
    //for (i = 0; i < row; ++i){
        //delete[] data[i];
        //delete[] dist[i];
    //}
    //delete[] data;
    //delete buffer;
    
    //}

    MPI_Finalize();
/*   
    // mds
    gsl_matrix* d = gsl_matrix_alloc(row, row);
    FILE* f = fopen("dist_matrix", "r");
    gsl_matrix_fscanf(f, d);
    fclose(f);
    MDS *mds = new MDS(d, 2);
    mds->scaling();
    gsl_matrix_free(d);
    delete mds;

*/

    return 0;
}
