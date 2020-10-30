/* File interpolation.h
 * Header for interpolation filledes */

#ifndef INTERPOLATION_H_ /* to ensure one single declaration - guard */
#define INTERPOLATION_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define R_OF_SUN    6.9634e10  /* cm */
#define ROWS        2482
typedef double ARRAY_POINTER[ROWS];  /* Add new variable type for arrays */
void   write_file(char file_name[], double x_arr[], int m, int n);
double kmtorsun(double x);
//void   read_file(void);
void read_file(ARRAY_POINTER RADIUS_FILE, double CSOUND_FILE[], 
     double DENSITY_FILE[], double PRESSURE_FILE[], double G1_FILE[],
     double TEMP_FILE[]);
void   linear_interp(double X[ROWS], double Y[ROWS], double a, double b, 
       int m, double X_i[], double Y_i[]);
//double RADIUS_FILE[ROWS], CSOUND_FILE[ROWS], DENSITY_FILE[ROWS], 
//       PRESSURE_FILE[ROWS], G1_FILE[ROWS], TEMP_FILE[ROWS];


#endif 
