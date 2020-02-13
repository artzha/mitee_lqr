//
//  matrix.h
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/9/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#ifndef matrix_h
#define matrix_h

#include <stdio.h>
#include <stdlib.h>

typedef struct MatrixStruct {
    int rows, cols;
    double **arr;
} Matrix;

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

/*
 Copies entries of arr into row or col of matrix, denoted by type parameter
 Setting type = 'R' fills row and type = 'C' fills col
 */
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

/*
 Multiplies matrix mult and matrix mcand together
 and returns product as new matrix
 
 Returns NULL on invalid matrix sizes
 */
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

/*
 Copies entries in matrix cpy into matrix mat beginning from start, end
 Returns:   1 for a successful copy and 0 for copying error
 */
int matrix_fill(Matrix *mat, Matrix*cpy, int start, int end) {
    int start_row = start;
    int start_col = end;
    for (int i = 0; i < cpy->rows; ++i) {
        if (start_row + i >= mat->rows) {
            return 0;
        } // check within row bounds
        for (int j = 0; j < cpy->cols; ++j) {
            if (start_col + j >= mat->cols) {
                return 0;
            } // check within col bounds
            mat->arr[start_row+i][start_col+j] = cpy->arr[i][j];
        }
    } // copies matrix entries from cpy into mat
    
    return 1;
}

#endif /* matrix_h */
