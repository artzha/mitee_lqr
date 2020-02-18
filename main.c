//
//  main.c
//  mitee-lqr
//
//  Created by Arthur Zhang on 2/9/20.
//  Copyright © 2020 Arthur Zhang. All rights reserved.
//

#include <stdio.h>
#include "matrix.h"

int main(int argc, const char * argv[]) {
	// testing inverse
	Matrix* mat = NULL;
	matrix_ctor(mat, 3, 3);

	mat->arr[0][0] = 1;
	mat->arr[0][1] = 2;
	mat->arr[0][2] = 3;

	mat->arr[1][0] = 2;
	mat->arr[1][1] = 2;
	mat->arr[1][2] = 2;

	mat->arr[2][0] = 4;
	mat->arr[2][1] = 5;
	mat->arr[2][2] = 1;

	Matrix* P = matrix_inverse(mat);

	matrix_dtor(mat);
	matrix_dtor(P);
    return 0;
}
