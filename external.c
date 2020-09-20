/*
 * external.c
 *
 *  Created on: Sep. 15, 2020
 *      Author: matblisc
 */

#include "external.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define NUM_SAMPLES 6867

static double data[NUM_SAMPLES][9];
static double torque[NUM_SAMPLES][3];


void loadSampleData() {
    FILE *file = fopen("sample-data.txt", "r");
    for (int i = 0; i < NUM_SAMPLES; i++) {
        fscanf(file,
               "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &(data[i][0]),
               &(data[i][1]),
               &(data[i][2]),
               &(data[i][3]),
               &(data[i][4]),
               &(data[i][5]),
               &(data[i][6]),
               &(data[i][7]),
               &(data[i][8]));
    } // for
    fclose(file);

    file = fopen("torque.txt", "r");
    for (int i = 0; i < NUM_SAMPLES; i++) {
        fscanf(file,
               "%lf %lf %lf",
               &(torque[i][0]),
               &(torque[i][1]),
               &(torque[i][2]));
    } // for
    fclose(file);
} // loadSampleData


void getAngularPosition(double* pos1, double* pos2, double* pos3) {
    static size_t i = 0;

    *pos1 = data[i][0];
    *pos2 = data[i][1];
    *pos3 = data[i][2];

    ++i;
} // getAngularPosition


void getAngularVelocity(double* vel1, double* vel2, double* vel3) {
    static size_t i = 0;

    *vel1 = data[i][3];
    *vel2 = data[i][4];
    *vel3 = data[i][5];

    ++i;
} // getAngularVelocity


void getMagneticField(double* b1, double* b2, double* b3) {
    static size_t i = 0;

    *b1 = data[i][6];
    *b2 = data[i][7];
    *b3 = data[i][8];

    ++i;
} // getMagneticField


void setMagnetorquer(double u1, double u2, double u3) {
    static size_t i = 0;

    printf("Timestep %ld\n", i);
    printf("Inputs: %e %e %e\n", u1, u2, u3);
    printf("Error:  %e %e %e\n\n",
            u1 + torque[i][0], // we should have torque = -u
            u2 + torque[i][1],
            u3 + torque[i][2]);

    ++i;
} // setMagnetorquer


// print matrix in human readable format
int printMatrix(FILE* stream, gsl_matrix* m, char* fmt) {
        size_t rows=m->size1;
        size_t cols=m->size2;
        size_t row,col,ml;
        int fill;
        char buf[100];
        gsl_vector *maxlen;

        maxlen=gsl_vector_alloc(cols);
        for (col=0;col<cols;++col) {
                ml=0;
                for (row=0;row<rows;++row) {
                        sprintf(buf,fmt,gsl_matrix_get(m,row,col));
                        if (strlen(buf)>ml)
                                ml=strlen(buf);
                }
                gsl_vector_set(maxlen,col,(double)ml);
        }

        for (row=0;row<rows;++row) {
                for (col=0;col<cols;++col) {
                        sprintf(buf,fmt,gsl_matrix_get(m,row,col));
                        fprintf(stream,"%s",buf);
                        fill=(int)gsl_vector_get(maxlen,col)+1-(int)strlen(buf);
                        while (--fill>=0)
                                fprintf(stream," ");
                }
                fprintf(stream,"\n");
        }
        gsl_vector_free(maxlen);
        fprintf(stream,"\n");
        return 0;
}
