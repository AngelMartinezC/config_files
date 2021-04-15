#include "pluto.h"
//#include "interpolation.h"

#define yres_ang    304         /* 4 pixels more than y resolution */
#define xres_ang    150         /* x resolution (not used) */
#define ROWS_ang    2482        /* Number of rows in solar model S */
#define COLS_ang    6           /* Number of columns S-model */
#define R_SUN       6.9634e10   /* Sun radius [cm] */
#define SOL_BOTTOM  1e4         /* Depth of simulation (as pluto_S.ini */
#define SOL_TOP     0.0         /* Photosphere */


void   linear_interp(double X[ROWS_ang], double Y[ROWS_ang], double a, 
       double b, int m, double X_i[], double Y_i[]);
int    find_index(double ARR[ROWS_ang], double value);
void   read_file(void);
double RADIUS_FILE[ROWS_ang], CSOUND_FILE[ROWS_ang], DENSITY_FILE[ROWS_ang],
       PRESSURE_FILE[ROWS_ang], G1_FILE[ROWS_ang], TEMP_FILE[ROWS_ang];
double derivative_b(double X[], double Y[], int idx);
double derivative_f(double X[], double Y[], int idx);
double derivative_c(double X[], double Y[], int idx);
void   hermite_interp(double X[ROWS_ang], double Y[ROWS_ang], double a, 
       double b, int m, double X_i[], double Y_i[]);
double psi_0(double z);
double psi_1(double z);



/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  v[PRS] = 1.0;
  v[RHO] = 10.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  v[TRC] = 0.0;
  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;
  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}



/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  int i, j, k;
  int id;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  int aa = 0;
  
  double s;
  read_file();
  double *RHO_i = malloc(yres_ang*sizeof(double));
  double *PRS_i = malloc(yres_ang*sizeof(double));
  double *R_i = malloc(yres_ang*sizeof(double));
  
  double UNIT_PRESS = UNIT_DENSITY*pow(UNIT_VELOCITY,2);
  double A = 1.0 - (SOL_BOTTOM*UNIT_LENGTH/R_SUN);  /* take care of units */
  double B = 1.0 - (SOL_TOP*UNIT_LENGTH/R_SUN);
  linear_interp(RADIUS_FILE, DENSITY_FILE,  A, B, yres_ang, R_i, RHO_i);
  linear_interp(RADIUS_FILE, PRESSURE_FILE,  A, B, yres_ang, R_i, PRS_i);
  TOT_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] = RHO_i[yres_ang-j]/UNIT_DENSITY;
    d->Vc[PRS][k][j][i] = PRS_i[yres_ang-j]/UNIT_PRESS;
  }
  free(R_i);
  free(RHO_i);
  free(PRS_i);

}




/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
{
}
#if PHYSICS == MHD



/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif



/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
}




#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}





/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
{
  return 0.0;
}
#endif




void read_file(void)
{
  int i;
  char ch;
  FILE *fp;
  double a, b, c, d, e, f;
  fp = fopen("Solar_Model.dat", "r");
  if (fp == NULL){
    perror("Error while reading file.\n");
    exit(EXIT_FAILURE);
  }
  while (fscanf(fp,"%lf %lf %lf %lf %lf %lf",&a,&b,&c,&d,&e,&f)==COLS_ang){
    RADIUS_FILE[i]    = a;
    CSOUND_FILE[i]    = b;
    DENSITY_FILE[i]   = c;
    PRESSURE_FILE[i]  = d;
    G1_FILE[i]        = e;
    TEMP_FILE[i]      = f;
    i++;
  }
  fclose(fp);
}



double derivative_b(double X[], double Y[], int idx){
  double up   = (Y[idx] - Y[idx-1]);
  double down = (X[idx] - X[idx-1]);
  return up/down;
}

double derivative_f(double X[], double Y[], int idx){
  double up   = (Y[idx+1] - Y[idx]);
  double down = (X[idx+1] - X[idx]);
  return up/down;
}

double derivative_c(double X[], double Y[], int idx){
  double up   = (Y[idx+1] - Y[idx-1]);
  double down = (X[idx+1] - X[idx-1]);
  return up/(2.0*down);
}


double psi_0(double z){
  double psi0 = 2*pow(z,3) - 3*pow(z,2) + 1;
  return psi0;
}
double psi_1(double z){
  double psi1 = pow(z,3) - 2*pow(z,2) + z;
  return psi1;
}


void   hermite_interp(double X[ROWS_ang], double Y[ROWS_ang], double a, 
       double b, int m, double X_i[], double Y_i[]){
  double dx = (b-a)/(double)m;
  double derivat[m];
  int idx;
  for (int k=0; k<m; k++){
    if (k==0){
      derivat[k] = derivative_f(X, Y, k);}
    else{
      derivat[k] = derivative_b(X, Y, k);}
    double s[10];
    for (int j=k; j<10; j+=1){
      s[j] = (X[k+1] - X[k])*j/10.0 + X[k];}
    double z[10];
    for (int i=0; i<10; i++){
      z[i] = (s[i] - X[k])/(X[k+1] - X[k]);}
    double *psi0, *psi1; 
    double apsi_0[10], apsi_1[10], H[10];
    for (int i=0; i<10; i++){
      H[i] = Y[k]*psi_0(z[i]) + Y[k+1]*psi_0(1-z[i]) +
             derivat[k]*X[k+1] - X[k]*psi_1(z[i]) - 
             derivat[k+1]*X[k+1] - X[k]*psi_1(1-z[i]);
    }

  }
}




/* Linear interpolation.
 */
void   linear_interp(double X[ROWS_ang], double Y[ROWS_ang], double a, 
       double b, int m, double X_i[], double Y_i[]){

  double dx = (b-a)/(double)m;
  double step,interp, error;
  double x0, x1, y0, y1, M, B;
  int idx0, idx1;
  FILE *fp;
  fp = fopen("Or.dat","w");
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
    //fprintf(fp,"%f");
    step   = step+dx;
  }
  fclose(fp);
}



/* Locate thenearest left-located index from value in array 
 */
int find_index(double ARR[ROWS_ang], double value){
  for (int i=0; i<ROWS_ang; i++){
    if (ARR[i]-value <= 0){
      return i;
    }
    else {}
  }
}
