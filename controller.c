//
//  controller.c
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/16/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#include "controller.h"

void initializeController(Controller *cntl) {
    /* Allocate all dynamic memory */
    /* Constants State Parameters */
    cntl->A_c = gsl_matrix_alloc(6, 6);
    cntl->A_d = gsl_matrix_alloc(6, 6);
    cntl->Q   = gsl_matrix_alloc(6, 6);
    cntl->R   = gsl_matrix_alloc(3, 3);
    cntl->J   = gsl_matrix_alloc(3, 3);
    /* Dynamic State Parameters */
    cntl->b    = gsl_vector_alloc(3);
    cntl->bmat = gsl_matrix_alloc(3, 3);
    cntl->B_c  = gsl_matrix_alloc(6, 3);
    cntl->B_d  = gsl_matrix_alloc(6, 3);
    cntl->P    = gsl_matrix_alloc(6, 6);
    cntl->K    = gsl_matrix_alloc(3, 6);
    /* Identity and Zero matrices */
    cntl->I_6    = gsl_matrix_alloc(6, 6);
    cntl->Zero_6 = gsl_matrix_alloc(6, 6);
    cntl->Zero_3 = gsl_matrix_alloc(3, 3);

    /* [kg*m^2] Approx for MiTEE 2 (From Three_MT main.m) */
    cntl->J12 = -0.7724;
    cntl->J31 = -0.0463;
    cntl->J23 = 0.7904;

    /* Initialize timestep */
    cntl->dt = 1; // 1 second

    /* TODO: Initialize equilibrium point */
    
    /* Initialize A_c matrix */
    double n                = cntl->eq;
    double neg_three_nn_J23 = -3*n*n*cntl->J23;
    double neg_n_J23        = -n*cntl->J23;
    double three_nn_J31     = 3*n*n*cntl->J31;
    double neg_n_J12        = -n*cntl->J12;

    double A_c = {0,                0,            n, 1,         0, 0,
                  0,                0,            0, 0,         1, 0,
                  -n,               0,            0, 0,         0, 1,
                  neg_three_nn_J23, 0,            0, 0,         0, neg_n_J23,
                  0,                three_nn_J31, 0, 0,         0, 0,
                  0,                0,            0, neg_n_J12, 0, 0};
    memcpy(cntl->A_c->data, A_c, 6*6*sizeof(double));

    /* Initialize A_d matrix */
    gsl_matrix_memcpy(cntl->A_d, cntl->A_c);
    gsl_matrix_scale(cntl->A_d, cntl->dt);
    gsl_linalg_exponential_ss(cntl->A_d, cntl->A_d, 0.01);

    /* Initialize J matrix */
    double J = {0.007, 0,      0,
                0,     0.0320, 0,
                0,     0,      0.0323};
    memcpy(cntl->J->data, J, 3*3*sizeof(double));

    /* Initialize identity and zero matrices */
    gsl_matrix_set_identity(cntl->I_6);
    gsl_matrix_set_zero(cntl->Zero_6);
    gsl_matrix_set_zero(cntl->Zero_3);

    /* Initialize Q (position) and R (inputs) cost matrices */
    double Q = {1.5e-7, 0,      0,      0,      0,      0,
                0,      1.5e-7, 0,      0,      0,      0,
                0,      0,      1.5e-7, 0,      0,      0,
                0,      0,      0,      1.5e-3, 0,      0,
                0,      0,      0,      0,      1.5e-3, 0,
                0,      0,      0,      0,      0,      1.5e-3};
    memcpy(cntl->Q->data, Q, 6*6*sizeof(double));
    
    double R = {1e7, 0,   0,
                0,   1e7, 0,
                0,   0,   1e7};
    memcpy(cntl->R->data, R, 6*6*sizeof(double));

    /* TODO: Initialize initial position */
    
} // initializeController

void computeDynamicInputs(Controller *cntl, int time) {
    // TODO: Figure out how to compute x position

    /* TODO: Poll magnetometers for magnetic field readings */

    /* Recompute B_d matrix */
    computeBMatrices(cntl);

    /* TODO: Update satellite position using sensor measurements */

    /* Compute P(t) matrix using recomputed matrices */
    computePMatrix(cntl);

    /* Compute K(t) gain matrix using P(t) */
    computeGainMatrix(cntl);

} // computeDynamicInputs

void computeBMatrices(Controller* cntl) {
    // declare all matrices
    static gsl_matrix* J_inv          = NULL;
    static gsl_matrix* quotient       = NULL;
    static gsl_matrix* J_inv_quotient = NULL;
    static gsl_matrix* A_c_inv        = NULL;
    static gsl_matrix* I_minus_A_d    = NULL;
    static gsl_matrix* transform      = NULL;

    // allocate memory on first function call
    if (!J_inv) {
        J_inv          = gsl_matrix_alloc(3, 3);
        quotient       = gsl_matrix_alloc(3, 3);
        J_inv_quotient = gsl_matrix_alloc(3, 3);
        A_c_inv        = gsl_matrix_alloc(6, 6);
        I_minus_A_d    = gsl_matrix_alloc(6, 6);
        transform      = gsl_matrix_alloc(6, 6);

        // compute matrices that are constant
        invert(cntl->J, J_inv);

        invert(cntl->A_c, A_c_inv);

        gsl_matrix_set_identity(I_minus_A_d);
        gsl_matrix_sub(I_minus_A_d, cntl->A_d);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, A_c_inv, I_minus_A_d, 0.0, transform);
    } // if
    
    /* Recompute B_c matrix */
    double bT_b;
    gsl_blas_ddot(cntl->b, cntl->b, bT_b);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1/bT_b, cntl->bmat, cntl->bmat, 0.0, quotient);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,    J_inv,      quotient,   0.0, J_inv_quotient);
    
    concatenate_vertical(cntl->Zero_3, J_inv_quotient, cntl->B_c);

    /* Recompute B_d matrix */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, transform, cntl->B_c, 0.0, cntl->B_d);
} // computeBdMatrix

void computePMatrix(Controller* cntl) {
    // P is computed using the process describe in Appendix B of the Sutherland paper
    // https://ieeexplore.ieee.org/document/8272473

    // declare all matrices
    static gsl_matrix* negativeQ        = NULL;
    static gsl_matrix* N                = NULL;
    static gsl_matrix* R_inverse        = NULL;
    static gsl_matrix* Rinv_BdT         = NULL;
    static gsl_matrix* Bd_Rinv_BdT      = NULL;
    static gsl_matrix* A_d_transpose    = NULL;
    static gsl_matrix* L                = NULL;
    static gsl_matrix* N_plus_L         = NULL;
    static gsl_matrix* N_plus_L_inverse = NULL;
    static gsl_matrix* N_minus_L        = NULL;
    static gsl_matrix* H                = NULL;
    static gsl_matrix* S                = NULL;
    static gsl_matrix* X1_inverse       = NULL;

    // allocate memory on first function call
    if (!negativeQ) {
        negativeQ        = gsl_matrix_alloc(6,  6);
        N                = gsl_matrix_alloc(12, 12);
        R_inverse        = gsl_matrix_alloc(3,  3);
        Rinv_BdT         = gsl_matrix_alloc(3,  6);
        Bd_Rinv_BdT      = gsl_matrix_alloc(6,  6);
        A_d_transpose    = gsl_matrix_alloc(6,  6);
        L                = gsl_matrix_alloc(12, 12);
        N_plus_L         = gsl_matrix_alloc(12, 12);
        N_plus_L_inverse = gsl_matrix_alloc(12, 12);
        N_minus_L        = gsl_matrix_alloc(12, 12);
        H                = gsl_matrix_alloc(12, 12);
        S                = gsl_matrix_alloc(12, 12);
        X1_inverse       = gsl_matrix_alloc(6,  6);

        // N matrix is constant
        gsl_matrix_memcpy(negativeQ, cntl->Q);
        gsl_matrix_scale(negativeQ, -1.0);
        concatenate2x2(cntl->A_d, cntl->Zero_6, negativeQ, cntl->I_6);

        // R^-1 and A_d^T are constant
        invert(cntl->R, R_inverse);
        gsl_matrix_transpose_memcpy(A_d_transpose, cntl->A_d);
    } // if

    // L matrix
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,   1.0, R_inverse, cntl->B_d, 0.0, Rinv_BdT);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cntl->B_d, Rinv_BdT,  0.0, Bd_Rinv_BdT);
    concatenate2x2(cntl->I_6, Bd_Rinv_BdT, cntl->Zero_6, A_d_transpose);

    // Hamiltonian matrix
    gsl_matrix_memcpy(N_plus_L, N);
    gsl_matrix_add(N_plus_L, L);
    invert(N_plus_L, N_plus_L_inverse);

    gsl_matrix_memcpy(N_minus_L, N);
    gsl_matrix_sub(N_minus_L, L);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, N_plus_L_inverse, N_minus_L, 0.0, H);

    // use Newton-Raphson process to determine S (positive sqrt of H)
    runNewtonRaphsonProcess(cntl, H, S);

    // extract X1 and X2 from H - S
    gsl_matrix_sub(H, S);
    gsl_matrix_view X1 = gsl_matrix_submatrix(H, 0, 0, 6, 6);
    gsl_matrix_view X2 = gsl_matrix_submatrix(H, 6, 0, 6, 6);

    // calculate P(t) = X2 * X1^-1
    invert(&X1.matrix, X1_inverse);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &X2.matrix, X1_inverse, 0.0, cntl->P);
} // computePMatrix

/* Helper function for computing P(t) */
void runNewtonRaphsonProcess(Controller* cntl, gsl_matrix* H, gsl_matrix* S) {
    // declare all matrices
    static gsl_matrix* S_prev    = NULL;
    static gsl_matrix* S_inverse = NULL;
    static gsl_matrix* H_squared = NULL;

    // allocate memory on first function call
    if (!S_prev) {
        S_prev          = gsl_matrix_alloc(12, 12);
        S_inverse       = gsl_matrix_alloc(12, 12);
        H_squared       = gsl_matrix_alloc(12, 12);

        // S should start as identity
        gsl_matrix_set_identity(S);
    } // if

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H, H, 0.0, H_squared);

    do {
        gsl_matrix_memcpy(S_prev, S);
        invert(S, S_inverse);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 0.5, S_inverse, H_squared, 0.5, S);
    } while(!newtonRaphsonConverged(S, S_prev));
} // runNewtonRaphsonProcess

// check if S matrix from Newton-Raphson iteration has converged
bool newtonRaphsonConverged(gsl_matrix* S, gsl_matrix* S_prev) {
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            if (gsl_matrix_get(S, i, j) - gsl_matrix_get(S_prev, i, j)) {
                return false;
            } // if
        } // for j
    } // for i

    return true;
} // newtonRaphsonConverged

/* Computes gain matrix K for magnetorquers*/
void computeGainMatrix(Controller* cntl) {
    // declare all matrices
    static gsl_matrix* P_Ad                = NULL;
    static gsl_matrix* BdT_P_Ad            = NULL;
    static gsl_matrix* P_Bd                = NULL;
    static gsl_matrix* BdT_P_Bd            = NULL;
    static gsl_matrix* R_plus_BdT_P_Bd     = NULL;
    static gsl_matrix* R_plus_BdT_P_Bd_inv = NULL;

    // allocate memory on first function call
    if (!P_Ad) {
        P_Ad                = gsl_matrix_alloc(6, 6);
        BdT_P_Ad            = gsl_matrix_alloc(3, 6);
        P_Bd                = gsl_matrix_alloc(6, 3);
        BdT_P_Bd            = gsl_matrix_alloc(3, 3);
        R_plus_BdT_P_Bd     = gsl_matrix_alloc(3, 3);
        R_plus_BdT_P_Bd_inv = gsl_matrix_alloc(3, 3);
    } // if

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cntl->P,   cntl->A_d, 0.0, P_Ad);
    gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, cntl->B_d, P_Ad,      0.0, BdT_P_Ad);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cntl->P,   cntl->B_d, 0.0, P_Bd);
    gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, cntl->B_d, P_Bd,      0.0, BdT_P_Bd);

    gsl_matrix_memcpy(R_plus_BdT_P_Bd, cntl->R);
    gsl_matrix_add(R_plus_BdT_P_Bd, BdT_P_Bd);
    invert(R_plus_BdT_P_Bd, R_plus_BdT_P_Bd_inv);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, R_plus_BdT_P_Bd_inv, BdT_P_Ad, 0.0, cntl->K);
} // computeGainMatrix

void sendMTInputs(/* Inputs */) {
    /* TODO: Compute new inputs to magnetorquers using gain matrix K */
    
    /* TODO: Send inputs to magnetorquers */
    
} // sendMTInputs

// concatenate four matrices into a single larger matrix
void concatenate2x2(gsl_matrix* a, gsl_matrix* b, gsl_matrix* c, gsl_matrix* d, gsl_matrix* result) {
    gsl_matrix_view aview = gsl_matrix_submatrix(result, 0,        0,        a->size1, a->size2);
    gsl_matrix_view bview = gsl_matrix_submatrix(result, 0,        a->size2, b->size1, b->size2);
    gsl_matrix_view cview = gsl_matrix_submatrix(result, a->size1, 0,        c->size1, c->size2);
    gsl_matrix_view dview = gsl_matrix_submatrix(result, a->size1, a->size2, d->size1, d->size2);

    gsl_matrix_memcpy( &aview.matrix, a);
    gsl_matrix_memcpy( &bview.matrix, b);
    gsl_matrix_memcpy( &cview.matrix, c);
    gsl_matrix_memcpy( &dview.matrix, d);
} // concatenate2x2

// concatenate two matrices vertically
void concatenate_vertical(gsl_matrix* a, gsl_matrix* b, gsl_matrix* result) {
    gsl_matrix_view aview = gsl_matrix_submatrix(result, 0,        0, a->size1, a->size2);
    gsl_matrix_view bview = gsl_matrix_submatrix(result, a->size1, 0, b->size1, b->size2);

    gsl_matrix_memcpy( &aview.matrix, a);
    gsl_matrix_memcpy( &bview.matrix, b);
} // concatenate2x2

// invert a matrix using LU decomposition
void invert(gsl_matrix* mat, gsl_matrix* inv) {
    // LU decomposition
    gsl_matrix* LU = gsl_matrix_alloc(mat->size1, mat->size2);
    gsl_matrix_memcpy(LU, mat);
    gsl_permutation* p = gsl_permutation_alloc(mat->size1);
    int signum = 0;
    gsl_linalg_LU_decomp(LU, p, &signum);

    // invert
    gsl_linalg_LU_invert(LU, p, inv);

    // deallocate memory
    gsl_matrix_free(LU);
    gsl_permutation_free(p);
} // invert

