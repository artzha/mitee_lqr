//
//  controller.c
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/16/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#include "controller.h"

void initializeController(Controller *cntl) {
    /* J constants taken from controller */
    
    /* Initialize all matrices to be used in detumbling procedure */
    matrix_ctor(&cntl->B_b, 3, 3);  // [~B_b] matrix
    matrix_ctor(&cntl->A_c, 6, 6);  // Ac matrix
    matrix_ctor(&cntl->B_c, 6, 3);  // Bc matrix
    matrix_ctor(&cntl->A_d, 6, 6);  // Ad matrix
    matrix_ctor(&cntl->B_d, 6, 3);  // Bd matrix
    matrix_ctor(&cntl->I_6 , 6, 6); // 6x6 identity matrix
    // TODO: Initialize identity matrix properly
        
    computeConstants(cntl);
    computeDynamicInputs(cntl, 0); // assume start at time 0
}

void computeConstants(Controller *cntl) {
    /* [kg*m^2] Approx for MiTEE 2 (From Three_MT main.m) */
    cntl->J12   = -0.7724;
    cntl->J13   = -0.0463;
    cntl->J23   = 0.7904;
    
    /* TODO: Ask about what to initilize Bb as */
    cntl->B_b_prev[0] = -3.211*0.000001;
    cntl->B_b_prev[1] = -2.378*0.00001;
    cntl->B_b_prev[2] = -2.026*0.000001;
        
    /* Discretize matrix A_c */
    matrix_scal_mult(&cntl->A_c, TIMESTEP);
    cntl->A_d = matrix_exp(&cntl->A_c);
}

void computeDynamicInputs(Controller *cntl, int time) {
    /* Discretize matrix B_c */
    
    // Step 1: Compute -A^-1
    Matrix *A_c_inv = matrix_inv(&cntl->A_c);
    matrix_scal_mult(A_c_inv, -1);
    // Step 2: Compute (I_6 - A_d)
    Matrix *I_sub_A_d = matrix_sub(&cntl->I_6, cntl->A_d);
    // Step 3: Compute B_c at time t
}
