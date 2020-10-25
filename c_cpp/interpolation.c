/* Interpolation I for data files
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ROWS     2482        /* Number of rows for Model S file */
#define COLS     6           /* Number of cols for Model S file */
#define PI       3.14159265  /* PI number because AJA */
#define R_SUN    6.9634e10   /* Sun radius [cm] */


double function(double x);
double linear(double a, double b, int m);
void   read_file(void);
int    find_index(double ARR[ROWS], double value);

double RADIUS[ROWS], CSOUND[ROWS], DENSITY[ROWS], PRESSURE[ROWS],
       G1[ROWS], TEMP[ROWS];


int main(char *argv[]){

  double y;
  //y = linear(3.9, 3.3, 11);
  read_file();
  printf("%f\n",RADIUS[0]);

  int val;

}




double function(double x){
  /* Dummy function for testing from CA */
  return 1/(25*pow(x,2)+1);
}



double linear(double a, double b, int m){
  /* Linear stepwise interpolation */
  double dx = (b-a)/(double)m;
  double x0, x1, M, B;
  double line;
  int val0, val1;
  for (int k=0; k<m; k++){
    val0 = find_index(RADIUS,0.5);
    val1 = find_index(RADIUS,1.0);
    if (val1-val0 == 0){
      // INTERPOLATE
    }
    else {
      // WRITE x0
    }
    //printf("k is %d\n",k);
    x0 = dx*k + a;
    x1 = dx*(k+1) + a;
    //M = (f(x0+dx)-f(x0))/(x1-x0);
    //B = f(x0) - M*x0;
    line = 2;
  }
}


int find_index(double ARR[ROWS], double value){

  for (int i=0; i<ROWS; i++){
    if (ARR[i]-value <= 0){
      return i;
    }
    else {}
  }
}


/* Read thermodynamical variables from the Solar Standard Model of 
 * Christensen Dalsgaard for a mangetic-free region.
 */

void read_file(void){

  int i;
  char ch;
  FILE *fp;
  double a, b, c, d, e, f;

  fp = fopen("Solar_Model.dat", "r");
  if (fp == NULL){
    perror("Error while reading file.\n");
    exit(EXIT_FAILURE);
  }

  while (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &a,&b,&c,&d,&e,&f) == COLS){
    RADIUS[i]     = a;
    CSOUND[i]    = b;
    DENSITY[i]   = c;
    PRESSURE[i]  = d;
    G1[i]        = e;
    TEMP[i]      = f;
    i++;
  }
  fclose(fp);
}




