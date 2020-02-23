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
    
    /* Initialize A_c matrix */
    
    /* Initialize A_d matrix */
    
    /* Initialize timestep */
    
    /* Initialize Q (position) and R (inputs) cost matrices */
    
    /* Initialize initial position */
    
}

void computeDynamicInputs(Controller *cntl, int time) {
    // TODO: Figure out how to compute x position
    
    /* Poll magnetometers for magnetic field readings */
    
    /* Recompute B_c matrix */
    
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

