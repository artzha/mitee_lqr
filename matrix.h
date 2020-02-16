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

/*
    TODO: Longterm Goals (in order of priority) :
        1. Restructure matrix ADT to represent sparse matrices
            more efficiently
        2. Optimize matrix by matrix operations
 */

typedef struct MatrixStruct {
    int rows, cols;
    double **arr;
} Matrix;

void matrix_ctor(Matrix* mat, int r, int c);

void matrix_dtor(Matrix *mat);

/*
 TODO: Multiplies all entries in matrix by scalar, saves result back
        into original matrix
 */
void matrix_scal_mult(Matrix *mat, double scalar);

/* Edits specific entry in matrix with value */
void matrix_fill(Matrix *mat, int row, int col, double value) {
    if (row >= mat->rows || col >= mat->cols) return;
    mat->arr[row][col] = value;
}

/*
 TODO: Raises matrix to the exponent and results result in output matrix
 */
Matrix* matrix_exp(Matrix *mat);

/*
 Copies entries of arr into row or col of matrix, denoted by type parameter
 Setting type = 'R' fills row and type = 'C' fills col
 */
int matrix_fill_rc(Matrix *mat, double *arr, int arr_sz, int start, char type);

/*
 Multiplies matrix mult and matrix mcand together
 and returns product as new matrix
 
 Returns NULL on invalid matrix sizes
 */
Matrix* matrix_mult(Matrix *mult, Matrix *mcand);

/*
 TODO: Adds two matrices
 */
Matrix* matrix_add(Matrix *lhs, Matrix*rhs);

/*
 TODO: Subtracts two matrices
 */
Matrix* matrix_sub(Matrix *lhs, Matrix*rhs);

/*
    TODO: Inverts matrix and returns result as new matrix
 */
Matrix* matrix_inv(Matrix *mat);

/*
 Copies entries in matrix cpy into matrix mat beginning from start, end
 Returns:   1 for a successful copy and 0 for copying error
 */
int matrix_cpy(Matrix *mat, Matrix*cpy, int start, int end);

#endif /* matrix_h */
