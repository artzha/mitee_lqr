/*
 * external.h
 *
 *  Created on: Sep. 15, 2020
 *      Author: matblisc
 *
 *  Functions for debugging that are used in place of interfacing
 *  with actual hardware during normal operation
 */

#ifndef EXTERNAL_H_
#define EXTERNAL_H_

#include <gsl/gsl_matrix.h>

// load sample data from MATLAB simulation
void loadSampleData();
void getAngularPosition(double *pos1, double *pos2, double *pos3);
void getAngularVelocity(double *vel1, double *vel2, double *vel3);
void getMagneticField(double *b1, double *b2, double *b3);

// print magnetorquer inputs
void setMagnetorquer(double u1, double u2, double u3);

// print GSL matrix
int printMatrix(FILE* stream, gsl_matrix* m, char* fmt);

#endif /* EXTERNAL_H_ */
