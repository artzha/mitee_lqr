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
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct ControllerState {
    /* Constants State Parameters */
    double eq; // equilibrium pt
    int dt;
    double J12, J23, J31;
    double NR_tolerance; // acceptable tolerance for Newton-Raphson process
    gsl_matrix* A_c;
    gsl_matrix* A_d;
    gsl_matrix* Q;
    gsl_matrix* R;
    gsl_matrix* J;
    /* Dynamic State Parameters */
    gsl_vector* b;
    gsl_matrix* bmat;
    gsl_matrix* B_c;
    gsl_matrix* B_d;
    gsl_matrix* P;
    gsl_matrix* K;
    /* Identity and Zero matrices */
    gsl_matrix* I_6;
    gsl_matrix* Zero_6;
    gsl_matrix* Zero_3;
} Controller;

typedef struct SensorState {
    double mag_x, mag_y, mag_z;
    double gyro_x, gyro_y, gyro_z;
} Sensors;

/* Declares all matrices and constants needed for LQR controller*/
void initializeController(Controller *cntl);

/* Computes values for numbers as discrete time steps */
void computeDynamicInputs(Controller *cntl, int time);

/* Computes B_c(t) and B_d(t) matrices */
void computeBMatrices(Controller* cntl);

/* Runs Newton Rhapson process to compute P(t) matrix*/
void computePMatrix(Controller* cntl);

/* Helper function for compute P(t) */
void runNewtonRaphsonProcess(Controller* cntl, gsl_matrix* H, gsl_matrix* S);

// check if S matrix from Newton-Raphson iteration has converged
bool newtonRaphsonConverged(gsl_matrix* S, gsl_matrix* S_prev);

/* Computes gain matrix K for magnetorquers*/
void computeGainMatrix(Controller* cntl);

/* Sends optimal inputs to magnetorquers for stabilization procedure */
void sendMTInputs(/* Inputs */);

// concatenate four matrices into a single larger matrix
void concatenate2x2(gsl_matrix* a, gsl_matrix* b, gsl_matrix* c, gsl_matrix* d, gsl_matrix* result);

// concatenate two matrices vertically
void concatenate_vertical(gsl_matrix* a, gsl_matrix* b, gsl_matrix* result);

// invert a matrix using LU decomposition
void invert(gsl_matrix* mat, gsl_matrix* inv);

#endif /* controller_h */
