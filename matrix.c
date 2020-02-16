//
//  matrix.c
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/9/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#include "matrix.h"

void matrix_ctor(Matrix* mat, int r, int c) {
    mat->rows = r;
    mat->cols = c;
    
    /* Allocate memory for matrix in dynamic memory */
    mat->arr = (double **) malloc(r * sizeof(double *));
    for (int i=0; i < r; i++) {
        mat->arr[i] = (double *)malloc(c * sizeof(double));
    }
}

void matrix_dtor(Matrix *mat) {
    /* Free dynamic memory for arr once finished */
    for (int i = 0; i < mat->rows; ++i) {
        free(mat->arr[i]);
    }
    free(mat->arr);
}

int matrix_fill_rc(Matrix *mat, double *arr, int arr_sz, int start, char type) {
    if (type == 'R' && arr_sz > mat->cols) return 0;
    if (type == 'C' && arr_sz > mat->rows) return 0;
    
    if (type == 'R') {
        for (int i = 0; i < arr_sz; ++i) {
            mat->arr[start][i] = arr[i];
        }
    } else if (type == 'C'){
        for (int i = 0; i < arr_sz; ++i) {
            mat->arr[i][start] = arr[i];
        }
    } else {
        return 0;
    }
    
    return 1;
}

Matrix* matrix_mult(Matrix *mult, Matrix *mcand) {
    if (mult->cols != mcand->rows) return NULL;
    
    /* Initialize return matrix*/
    Matrix* mat = NULL;
    matrix_ctor(mat, mult->rows, mcand->cols);
    
    /* Compute matrix product*/
    for (int i = 0; i < mult->rows; ++i) {
        for (int j = 0; j < mult->cols; ++j) {
            double running_sum = 0;
            for (int k = 0; k < mat->rows; ++k) {
                running_sum += mult->arr[i][k]*mcand->arr[k][i];
            } // compute sum of row i col j
            mat->arr[i][j] = running_sum;
        } // loop through each column
    } // loop through each row
    
    return mat;
}




