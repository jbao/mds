// $Id: pca.c 246 2011-05-09 14:32:05Z jbao $

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

gsl_matrix* diag_alloc(gsl_vector* X){
    gsl_matrix* mat = gsl_matrix_calloc(X->size, X->size);
    gsl_vector_view diag = gsl_matrix_diagonal(mat);
    gsl_vector_memcpy(&diag.vector, X);
    return mat;
}

void print2file(char* filename, gsl_matrix* mat){
    FILE* f = fopen(filename, "w");
    int i, j;
    for (i = 0; i < mat->size1; ++i){
        for (j = 0; j < mat->size2; ++j){
            fprintf(f, "%f ", gsl_matrix_get(mat, i, j));
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

int main(int argc, char **argv){
    int row = atoi(argv[2]);
    int col = atoi(argv[3]);
    printf("%d %d\n", row, col);
    gsl_matrix* data = gsl_matrix_alloc(row, col);
    //gsl_matrix* data = gsl_matrix_alloc(col, row);
    FILE* f = fopen(argv[1], "r");
    gsl_matrix_fscanf(f, data);
    //gsl_matrix_fread(f, data);
    //gsl_matrix_transpose_memcpy(data, data_raw);
    fclose(f);

    //printf("%f %f", gsl_matrix_get(data,0,0), gsl_matrix_get(data,0,1));
    //f = fopen("test.dat", "w");
    //gsl_matrix_fprintf(f, data, "%f");
    //fclose(f);
    
    // data centering, subtract the mean in each dimension (col.-wise)
    int i, j;
    double mean, sum, std;
    gsl_vector_view col_vector;
    for (i = 0; i < col; ++i){
        col_vector = gsl_matrix_column(data, i);
        mean = gsl_stats_mean((&col_vector.vector)->data, 1, 
            (&col_vector.vector)->size);
        gsl_vector_add_constant(&col_vector.vector, -mean);
        gsl_matrix_set_col(data, i, &col_vector.vector);
    }
    
    char filename[50];
    //sprintf(filename, "%s.zscore", argv[1]);
    //print2file(filename, data);

    gsl_matrix* u;
    if (col > row) {
        u = gsl_matrix_alloc(data->size2, data->size1);
        gsl_matrix_transpose_memcpy(u, data);
    } 
    else {
        u = gsl_matrix_alloc(data->size1, data->size2);
        gsl_matrix_memcpy(u, data);
    }

    // svd
    gsl_matrix* X = gsl_matrix_alloc(col, col);
    gsl_matrix* V = gsl_matrix_alloc(u->size2, u->size2);
    gsl_vector* S = gsl_vector_alloc(u->size2);
    gsl_vector* work = gsl_vector_alloc(u->size2);
    gsl_linalg_SV_decomp(u, V, S, work);
    //gsl_linalg_SV_decomp_jacobi(u, V, S);

    // mode coefficient
    //print2file("u.dat", u);
    /*
    // characteristic mode
    gsl_matrix* diag = diag_alloc(S);
    gsl_matrix* mode = gsl_matrix_alloc(diag->size1, V->size1);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diag, V, 0.0, mode);
    gsl_matrix_transpose(mode);
    print2file("mode.dat", mode);
    gsl_matrix_transpose(mode);
    */
    // reconstruction
    gsl_matrix *recons = gsl_matrix_alloc(u->size2, data->size1);
    if (col > row) {
        gsl_matrix_view data_sub = gsl_matrix_submatrix(data, 0, 0, 
            u->size2, u->size2);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, 
            &data_sub.matrix, 0.0, recons);
    }
    else
        gsl_blas_dgemm(CblasTrans, CblasTrans, 1.0, V, data, 0.0, recons);

    //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, u, mode, 0.0, 
    //    recons);
    gsl_matrix *recons_trans = gsl_matrix_alloc(recons->size2, 
        recons->size1);
    gsl_matrix_transpose_memcpy(recons_trans, recons);
    // take the first two eigenvectors
    gsl_matrix_view final = gsl_matrix_submatrix(recons_trans, 0, 0, 
            recons_trans->size1, 2);

    print2file(argv[4], &final.matrix);

    // eigenvalue
    gsl_vector_mul(S, S);
    f = fopen("eigenvalue.dat", "w");
    //gsl_vector_fprintf(f, S, "%f");
    fclose(f);

    gsl_matrix_free(data);
    gsl_matrix_free(X);
    gsl_matrix_free(V);
    //gsl_matrix_free(diag);
    //gsl_matrix_free(mode);
    gsl_matrix_free(recons);
    gsl_matrix_free(recons_trans);
    gsl_matrix_free(u);
    gsl_vector_free(S);
    gsl_vector_free(work);
    //gsl_vector_free(zero);
    //gsl_vector_free(corrcoef);
    //gsl_vector_free(corrcoef_mean);
    return 0;
}
