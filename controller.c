//
//  controller.c
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/16/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#include "controller.h"

void initializeController(Controller *cntl) {
    /* [kg*m^2] Approx for MiTEE 2 (From Three_MT main.m) */
    cntl->J12   = -0.7724;
    cntl->J13   = -0.0463;
    cntl->J23   = 0.7904;

    /* Initialize timestep */
    cntl->dt    =   1; // 1 second

    /* Initialize equilibrium point */
    
    /* Initialize A_c matrix */
    double n                = cntl->eq;
    double neg_three_nn_J23 =   -3*n*n*cntl->J23;
    double neg_n_J23        =   -n*cntl->J23;
    double three_nn_J31     =   3*n*n*cntl->J31;
    double neg_n_J12        =   -n*cntl->J12;

    double A_c      =  {0, 0, n, 1, 0, 0,
                        0, 0, 0,  0, 1, 0,
                        -n, 0, 0, 0, 0, 1,
                        neg_three_nn_J23, 0, 0, 0, 0, neg_n_J23,
                        0, three_nn_J31, 0, 0, 0, 0,
                        0, 0, 0, neg_n_J12, 0, 0};

    int col_sz      =   6;
    int row_sz      =   6;
    cntl->A_c       =   gsl_matrix_alloc(row_sz, col_sz);
    memcpy(cntl->A_c->data, A_c, col_sz*row_sz*sizeof(double));
    /* Initialize A_d matrix */
    memcpy(cntl->A_d->data, A_c, col_sz*row_sz*sizeof(double));
    gsl_matrix_scale(cntl->A_d, cntl->dt);
    gsl_linalg_exponential_ss(cntl->A_d, cntl->A_d, 0.01);
    
    /* Initialize top half of B_c matrix to 0s */
    double zeros_three =   {0, 0, 0};
    gsl_matrix_set_row(cntl->B_c, zeros_three, 0);
    gsl_matrix_set_row(cntl->B_c, zeros_three, 1);
    gsl_matrix_set_row(cntl->B_c, zeros_three, 2);

    /* Initialize J matrix */
    double J        =   {0.007, 0, 0,
                         0, 0.0320, 0,
                         0, 0, 0.0323};
    memcpy(cntl->J->data, J, col_sz*row_sz*sizeof(double));

    /* Initialize Q (position) and R (inputs) cost matrices */
    
    /* Initialize initial position */
    
}

void computeDynamicInputs(Controller *cntl, int time) {
    // TODO: Figure out how to compute x position
    
    /* Poll magnetometers for magnetic field readings */
    
    /* Recompute B_c matrix */
    gsl_matrix* J_inv;
    
    /* Recompute B_d matrix */
    
    /* Update satellite position using sensor measurements */
    
    /* Compute P(t) matrix using recomputed matrices */
    computePMatrix();
    
    /* Compute K(t) gain matrix using P(t) */
    
}

void computePMatrix(/* Inputs */) {
    // Calls helper function for numerically computing P(t)
}

void sendMTInputs(/* Inputs */) {
    /* Compute new inputs to magnetorquers using gain matrix K */
    
    /* Send inputs to magnetorquers */
    
}

