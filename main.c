/* main.c
 *
 *  Created on: Feb 20, 2020
 *      Author: matblisc
 */

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

#include "controller.h"

int main() {
    int n = 3;

    double A_data[] = {1, 2, 3,
                       0, 1, 4,
                       5, 6, 0};

    double B_data[] = {3, 0, 0,
                       0, 3, 0,
                       0, 0, 3};

    gsl_matrix* A = gsl_matrix_alloc(n, n);
    gsl_matrix* B = gsl_matrix_alloc(n, n);
    memcpy(A->data, A_data, n*n*sizeof(double));
    memcpy(B->data, B_data, n*n*sizeof(double));

    gsl_matrix* LU = gsl_matrix_alloc(n, n);
    gsl_matrix_memcpy(LU, A);
    gsl_permutation* p = gsl_permutation_alloc(n);
    int signum = 0;
    gsl_linalg_LU_decomp(LU, p, &signum);

    gsl_matrix* inverse = gsl_matrix_alloc(n, n);
    gsl_linalg_LU_invert(LU, p, inverse);

    gsl_matrix* product = gsl_matrix_alloc(n, n);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, inverse, 0.0, product);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%f ", gsl_matrix_get(product, i, j));
        }
        printf("\n");
    }


    return 0;
}
