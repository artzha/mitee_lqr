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

// timestep to used for descretized detumbling procedure
#define TIMESTEP    4

typedef struct ControllerState {
    int dt;
    double J12, J23, J13;
    double B_b_prev[3];
//    Matrix B_b;
//    Matrix I_6; // 6x6 identity matrix
//    Matrix A_c, B_c;
//    Matrix *A_d, *B_d;
    
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
