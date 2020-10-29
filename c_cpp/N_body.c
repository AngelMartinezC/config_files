/* Runge Kutta 4 method for solving a 2nd order ODE */

#include <stdio.h> 
#include <math.h>   
#include <stdlib.h> 

#define PI        3.141592654
#define G         6.67e-11       /* m^3 kg^-1 s^-1 */
#define MSUN      1.989e30       /* kg */
#define MEARTH    5.972e24       /* kg */
#define MJUP      1.898e27       /* kg */
#define MMARS     6.418e23       /* kg */
#define MSAT      5.683e26       /* kg */
#define MURAN     8.686e25       /* kg */
#define MNEPT     1.024e26       /* kg */
#define NPART     7              /* Number of particles */

double f1(double t, double X[], double Y[], double VX[], double VY[],
       double M[], double D1X[], double D1Y[]);
double f2(double t, double X[], double Y[], double VS[], double VY[],
       double M[], double D2X[], double D2Y[]);
void rk4_2(double(*d1x)(double,double*,double*,double*,double*,double*,
       double*,double*), double(*d2x)(double,double*,double*,double*,
       double*,double*,double*,double*), double ti, double X[], double Y[], 
       double V[], double U[], double M[], double tf, double Xf[], 
       double Yf[], double Vf[], double Uf[]);

/* int main: donde se hace el codigo en c++ */
int main(){

  double ti, tf, dt, tmax;
  double xi, xf, yi, yf, vi, vf, ui, uf;
  FILE *fp;

  double Xf[NPART], Yf[NPART], Vf[NPART], Uf[NPART];
  double M[NPART] = {MEARTH, MSUN, MJUP, MMARS, MSAT, MURAN, MNEPT};
  double X[NPART] = {1.5e11, 0.0, 7.66e11, 2.27e11, 14.33e11, 28.5e11, 
                     45.0e11};
  double Y[NPART] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double V[NPART] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double U[NPART] = {2.9e4, 0.0, 1.3e4, 2.4e4, 9.6e3, 6.8e3, 5.43e3};
  ti = 0.0;        /* Tiempo inicial de integration */             
  dt = 24.0*10*3600.0;     /* Temporal resolution */
  tmax = 42.0*365.25*24.0*3600.0;     /* tiempo final de integration */
 
  int aa = 0, max_step= (int) (tmax/dt);
  fp = fopen("DATA.dat","w");
  while (ti <= tmax){
   tf = ti + dt; 
   printf("Step number %d of %d\n",aa, max_step);
   rk4_2(f1,f2,ti,X,Y,V,U,M,tf,Xf,Yf,Vf,Uf);  /* Hacer RK4 */
   fprintf(fp,"%f ", tf);
   for (int i=0; i<NPART; i++){
     fprintf(fp,"%f %f ", Xf[i],Yf[i]);
   }
   fprintf(fp,"\n");
   ti = tf;   /* Actualizacion para paar al siguiente paso temporal */
   for (int i=0; i<NPART; i++){
     X[i] = Xf[i];     /* Actualizacion de la posicion */
     Y[i] = Yf[i];     /* Actualizacion de la posicion */
     V[i] = Vf[i];     /* Actualizacion de la posicion */
     U[i] = Uf[i];     /* Actualizacion de la posicion */
   }
   aa += 1;
  }
  fclose(fp);
}


/* Function 1: La que tiene en cuenta la velocidad */
double f1(double t, double X[], double Y[], double VX[], double VY[], 
       double M[], double D1X[], double D1Y[]){
  for (int i=0; i<NPART; i++){
    D1X[i] = VX[i];
    D1Y[i] = VY[i];
  }
}

double f2(double t, double X[], double Y[], double VX[], double VY[],
       double M[], double D2X[], double D2Y[]){
  long double pot = 0.0;
  for (int i=0; i<NPART; i++){
    double XX=0, YY=0;
    for (int j=0; j<NPART; j++){
      if (j!=i){
        pot = G*M[i]*M[j]/pow(pow(X[j]-X[i],2)+pow(Y[j]-Y[i],2),3.0/2.0);
        XX += pot*(X[j]-X[i]);
        YY += pot*(Y[j]-Y[i]);
      }
    }
    D2X[i] = (1.0/M[i])*XX;
    D2Y[i] = (1.0/M[i])*YY;
  }
}


/* RK4 solver for a 2nd order ODE */
void rk4_2(double(*d1x)(double,double*,double*,double*,double*,double*,
       double*,double*), double(*d2x)(double,double*,double*,double*,
       double*,double*,double*,double*), double ti, double X[], double Y[], 
       double V[], double U[], double M[], double tf, double Xf[], 
       double Yf[], double Vf[], double Uf[]){

  double h, t;
  double xxi[NPART], yyi[NPART], vvi[NPART], uui[NPART];
  double *K1X=malloc(NPART*sizeof(double)),*K1Y=malloc(NPART*sizeof(double));
  double *K1V=malloc(NPART*sizeof(double)),*K1U=malloc(NPART*sizeof(double));
  double *K2X=malloc(NPART*sizeof(double)),*K2Y=malloc(NPART*sizeof(double));
  double *K2V=malloc(NPART*sizeof(double)),*K2U=malloc(NPART*sizeof(double));
  double *K3X=malloc(NPART*sizeof(double)),*K3Y=malloc(NPART*sizeof(double));
  double *K3V=malloc(NPART*sizeof(double)),*K3U=malloc(NPART*sizeof(double));
  double *K4X=malloc(NPART*sizeof(double)),*K4Y=malloc(NPART*sizeof(double));
  double *K4V=malloc(NPART*sizeof(double)),*K4U=malloc(NPART*sizeof(double));
  h =tf-ti;  /* Rango temporal de integration */
  t=ti;
  d1x(t, X, Y, V, U, M, K1X, K1Y);
  d2x(t, X, Y, V, U, M, K1V, K1U);

  for (int i=0; i<NPART; i++){
    xxi[i] = X[i] + h*K1X[i]/2.0;
    yyi[i] = Y[i] + h*K1Y[i]/2.0;
    vvi[i] = V[i] + h*K1V[i]/2.0;
    uui[i] = U[i] + h*K1U[i]/2.0;
    //printf("Xi %f  Yi %f  ---- i %d\n",X[i],Y[i],i);
  }
  d1x(t+h/2.0, xxi, yyi, vvi, uui, M, K2X, K2Y);
  d2x(t+h/2.0, xxi, yyi, vvi, uui, M, K2V, K2U);

  for (int i=0; i<NPART; i++){
    xxi[i] = X[i] + h*K2X[i]/2.0;
    yyi[i] = Y[i] + h*K2Y[i]/2.0;
    vvi[i] = V[i] + h*K2V[i]/2.0;
    uui[i] = U[i] + h*K2U[i]/2.0;
  }
  d1x(t+h/2.0, xxi, yyi, vvi, uui, M, K3X, K3Y);
  d2x(t+h/2.0, xxi, yyi, vvi, uui, M, K3V, K3U);

  for (int i=0; i<NPART; i++){
    xxi[i] = X[i] + h*K3X[i];
    yyi[i] = Y[i] + h*K3Y[i];
    vvi[i] = V[i] + h*K3V[i];
    uui[i] = U[i] + h*K3U[i];
  }
  d1x(t+h, xxi, yyi, vvi, uui, M, K4X, K4Y);
  d2x(t+h, xxi, yyi, vvi, uui, M, K4V, K4U);

  for (int i=0; i<NPART; i++){
    Xf[i] = X[i] + (h*K1X[i] + 2.0*h*(K2X[i]+K3X[i]) + h*K4X[i])/6.0;
    Yf[i] = Y[i] + (h*K1Y[i] + 2.0*h*(K2Y[i]+K3Y[i]) + h*K4Y[i])/6.0;
    Vf[i] = V[i] + (h*K1V[i] + 2.0*h*(K2V[i]+K3V[i]) + h*K4V[i])/6.0;
    Uf[i] = U[i] + (h*K1U[i] + 2.0*h*(K2U[i]+K3U[i]) + h*K4U[i])/6.0;
  }
}












