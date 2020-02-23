//
//  controller.h
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/16/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#ifndef controller_h
#define controller_h

#include <stdio.h>
#include <gsl/gsl_matrix.h>
<<<<<<< HEAD
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
=======

typedef struct ControllerState {
    /* Constants State Paramters */
    double eq; // equilibrium pt
    int dt;
    double J12, J23, J13;
    gsl_matrix* A_c;
    gsl_matrix* A_d;
    gsl_matrix* Q;
    gsl_matrix* R;
    gsl_matrix* J;
    /* Dynamic State Parameters */
    gsl_matrix* B_b;
    gsl_matrix* B_c;
    gsl_matrix* B_d;
    gsl_matrix* P;
    gsl_matrix* K;
    gsl_matrix* I_6, Zero_6; // 6x6 identity/zero matrices
} Controller;

typedef struct SensorState {
    double mag_x, mag_y, mag_z;
    double gyro_x, gyro_y, gyro_z;
} Sensors;

/* Declares all matrices and constants needed for LQR controller*/
void initializeController(Controller *cntl);

/* Computes values for numbers as discrete time steps */
void computeDynamicInputs(Controller *cntl, int time);

/* Runs Newton Rhapson process to compute P(t) matrix*/
void computePMatrix(/* Inputs */);

/* Helper function for compute P(t) */
void runNewtonRhapsonProcess(/* Inputs*/);

/* Computes gain matrix K for magnetorquers*/
void computeGainMatrix(/* Inputs */);

/* Sends optimal inputs to magnetorquers for stabilization procedure */
void sendMTInputs(/* Inputs */);

#endif /* controller_h */
