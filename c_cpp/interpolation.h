/* File interpolation.h
 * Header for interpolation filledes */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define R_OF_SUN    6.9634e10  /* cm */

extern void   write_file(char file_name[], double x_arr[], int m, int n);
extern double kmtorsun(double x);
int hola(int x);
