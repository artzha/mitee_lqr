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

// returns a new zero matrix of size n x m
Matrix* matrix_zero(int n, int m) {
    // initialize return matrix
    Matrix* mat = NULL;
    matrix_ctor(mat, n, m);

    // fill all entries with zero
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            mat->arr[i][j] = 0;
        } // for j
    } // for i

    return mat;
} // matrix_zero

// returns a new identity matrix of size n
Matrix* matrix_identity(int n) {
    // start by initializing a zero matrix
    Matrix* mat = matrix_zero(n, n);

    // set the diagonals to 1
    for (int i = 0; i < n; ++i) {
        mat->arr[i][i] = 1;
    } // for i

    return mat;
} // matrix_identity

// adds matrices A and B together
// returns sum as a new matrix
// returns NULL on invalid matrix sizes
Matrix* matrix_add(Matrix* A, Matrix* B) {
    if (A->rows != B->rows || A->cols != B->cols) return NULL;

    // initialize return matrix
    Matrix* mat = NULL;
    matrix_ctor(mat, A->rows, A->cols);

    // compute sum
    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->rows; ++j) {
                mat->arr[i][j] = A->arr[i][j] + B->arr[i][j];
        } // for j
    } // for i

    return mat;
} // matrix_add

// subtracts matrix B from matrix A
// returns difference as a new matrix
// returns NULL on invalid matrix sizes
Matrix* matrix_sub(Matrix* A, Matrix* B) {
    if (A->rows != B->rows || A->cols != B->cols) return NULL;

    // initialize return matrix
    Matrix* mat = NULL;
    matrix_ctor(mat, A->rows, A->cols);

    // compute difference
    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->rows; ++j) {
            mat->arr[i][j] = A->arr[i][j] - B->arr[i][j];
        } // for j
    } // for i

    return mat;
} // matrix_sub

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

// computes the inverse of matrix A using LUP decomposition
// returns the inverse as a new matrix
// returns NULL if the matrix is non-invertible
Matrix* matrix_inverse(Matrix* A) {
    if(A->rows != A->cols) return NULL;

    int n = A->rows;

    // determine pivot matrix P
    Matrix* P = matrix_identity(A->rows);
    for (int i = 0; i < n; ++i) {
        // we need to find the largest (absolute value) element of each column at or below the diag
        double max_val = 0;
        int max_row = 0;
        for (int j = i; j < n; ++j) {
            // get the absolute value of the element
            int val = A->arr[j][i];
            if (val < 0) val *= -1;

            if (val > max_val) {
                max_val = val;
                max_row = i;
            } // if
        } // for j

        // if the rest of the column was zero, matrix is non-invertible
        if (max_val == 0) return NULL;

        // swap rows of P so that the diags of PA have the max elements
        if (max_row != i) {
            for (int j = 0; j < n; ++j) { // j is the column
                double temp = P->arr[i][j];
                P->arr[i][j] = P->arr[max_row][j];
                P->arr[max_row][j] = temp;
            } // for j
        } // if
    } // for i

    // Next step: find the LU decomposition of PA

    // initialize the matrix to store the result
    Matrix* LU = NULL;
    matrix_ctor(LU, n, n);

    // both L and U are stored in the same matrix mat since the only place where
    // both are nonzero is the diag and we know the diag of L is all ones

    for (int i = 0; i < n; ++i) {
        // calculate the Uij across row i
        for (int j = i; j < n; ++j) {
            // Uij = Aij - (product of U's above and L's to the left)
            LU->arr[i][j] = A->arr[i][j];
            for (int k = 0; k < i; ++k) {
                LU->arr[i][j] -= LU->arr[i][k] * LU->arr[k][j];
            } // for k
        } // for j

        // calculate the Lji down column i
        for (int j = i + 1; j < n; ++j) {
            // Lji = [Aji - (product of U's above and L's to the left)] / Uii
            LU->arr[j][i] = A->arr[j][i];
            for (int k = 0; k < i; ++k) {
                LU->arr[j][i] -= LU->arr[j][k] * LU->arr[k][i];
            } // for k
            LU->arr[j][i] /= LU->arr[i][i];
        } // for j
    } // for i

    // we must find the matrix X such that AX = I, or equivalently PAX = LUX = p
    // we will do this by solving PLUx = p where each column x of X and p of P
    // then we can simply solve Ly = p, and then use that solution to solve Ux = y
    Matrix* Y = NULL;
    matrix_ctor(Y, n, n);

    for (int i = 0; i < n; ++i) {
        // compute solution of Ly = p using forward substitution down the rows
        for (int j = 0; j < n; ++j) {
            // yji = pji - (sum of Ljk*yki)
            Y->arr[j][i] = P->arr[j][i];
            for (int k = 0; k < j; ++k) {
                Y->arr[j][i] -= LU->arr[j][k] * Y->arr[k][i];
            } // for k
        } // for j

        // compute solution of Ux = y using back substitution up the rows
        for (int j = n - 1; j >= 0; --j) {
            // xji = [yji - (sum of Ujk*xki)] / Ujj
            // we're done with P now, so we'll re-use it for X
            P->arr[j][i] = Y->arr[j][i];
            for (int k = 1; k < n; ++k) {
                P->arr[j][i] -= LU->arr[j][k] * P->arr[k][i];
            } // for k
            P->arr[j][i] /= LU->arr[j][j];
        } // for j
    } // for i

    matrix_dtor(LU);
    matrix_dtor(Y);
    return P;
} // matrix_inverse
