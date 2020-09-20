/* main.c
 *
 *  Created on: Feb 20, 2020
 *      Author: matblisc
 */

#include "controller.h"
#include "external.h"

#define NUM_LOOPS 6867

int main() {
    loadSampleData();

    Controller cntl;
    initializeController(&cntl);

    for (int i = 0; i < NUM_LOOPS; ++i) {
        /* Update measurements from attitude determination and magnetomers */
        updateSensors(&cntl);

        /* Recompute B_c and B_d matrices */
        computeBMatrices(&cntl);

        /* Compute P(t) matrix using recomputed matrices */
        computePMatrix(&cntl);

        /* Compute K(t) gain matrix using P(t) */
        computeGainMatrix(&cntl);

        /* Calculate and send magnetorquer inputs */
        sendMTInputs(&cntl);

    } // while

    return 0;
} // main
