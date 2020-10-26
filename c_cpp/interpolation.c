/* Interpolation I for data files
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ROWS     2482        /* Number of rows for Model S file */
#define COLS     6           /* Number of cols for Model S file */
#define PI       3.14159265  /* PI number because AJA */
#define R_SUN    6.9634e10   /* Sun radius [cm] */

typedef double ARRAY_POINT[ROWS];  /* Add new variable type for arrays */

double linear_interp(double X[ROWS], double Y[ROWS], double a, double b, 
       int m, double X_i[], double Y_i[]);
void   read_file(void);
int    find_index(double ARR[ROWS], double value);
void   write_file(char file_name[], double x_arr[], int m);

double RADIUS[ROWS], CSOUND[ROWS], DENSITY[ROWS], PRESSURE[ROWS],
       G1[ROWS], TEMP[ROWS];



int main(char *argv[]){

  double y;
  read_file();
  
  int m = 40;
  //ARRAY_POINT R_i, RHO_i;  /* When using pointers */
  double *RHO_i = malloc(m*sizeof(double));
  double *R_i = malloc(m*sizeof(double));

  //for (int i=0; i<=10; i++)
  printf("Pointer %f.  No pointer %f\n",RADIUS[0]);
  y = linear_interp(RADIUS, DENSITY, 0.9999, 1.0, m, R_i, RHO_i);

  printf("Pointer %f.  No pointer %f\n",R_i[0],RHO_i[m]*1.0e7);

  write_file("NEW.dbl",R_i, m);

}




/* Linear interpolation.
 */
double linear_interp(double X[ROWS], double Y[ROWS], double a, double b, 
       int m, double X_i[], double Y_i[]){

  double dx = (b-a)/(double)m;
  double step,interp, error;
  double x0, x1, y0, y1, M, B;
  int idx0, idx1;
  /* To ensure it is a monotonic crecent index */
  if (X[find_index(X,a)] < X[find_index(X,b)]){
    step = b;
    dx   = -1.0*dx;}
  else{
    step = a;}
  for (int k=0; k<=m; k++){
    idx0 = find_index(X,step);
    idx1 = idx0 + 1;
    x0  = X[idx0];
    x1  = X[idx1];
    y0  = Y[idx0];
    y1  = Y[idx1];
    M   = (y1-y0)/(x1-x0);
    B   = y0 - M*x0;
    interp  = M*step + B;
    error   = (Y[idx0]-interp)/Y[idx0];
    //(*X_i)[k] = step;
    //(*Y_i)[k] = step;
    X_i[k] = step;
    Y_i[k] = interp;
    step   = step+dx;
  }
}



/* Locate thenearest left-located index from value in array 
 */

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


/* Save data into file_name
 */

void write_file(char file_name[], double x_arr[], int m){

  FILE *fp;

  fp = fopen(file_name,"w");
  //for (int i=m; i>=0; i--){
  for (int i=0; i<=m; i++){
    printf("%f\n",x_arr[i]);
    fprintf(fp, "%f\n", x_arr[i]);
    }
  fclose(fp);
}


