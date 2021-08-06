/* Attempt to evolve solar interior conditions from the PLUTO module. 
 *
 * if change DEPTH, chande 100/6 for photosphere calibration (zero depth).
 *
 * Major changes:
 *  - Change endianity for InputDataOpen() MACRO --> little endianess
 *  - Read gravity data on runtime
 *  - P_TOT as lateral condition
 *
 */

#include "pluto.h"

#define  x1_left      g_domBeg[0]
#define  x1_right     g_domEnd[0]
#define  x2_top       g_domBeg[1]
#define  x2_bottom    g_domEnd[1]

double UNIT_PRESSURE = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;

double interp_rho(double y);
double interp_prs(double y);
double interp_gra(double y);
double rho_min;  /* Values taken from the solar model at depth 0 and 10 Mm */
double rho_max;
double prs_min;
double prs_max;


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  int first_call = 1;

  /* Call of upper and lower density, pressure variables */
  if (first_call){
    rho_min = interp_rho(x2_bottom);
    rho_max = interp_rho(x2_top);
    prs_min = interp_prs(x2_bottom);
    prs_max = interp_prs(x2_top);
    first_call = 0;
  }

  v[RHO] = interp_rho(x2);
  v[PRS] = interp_prs(x2);
  v[VX1] = v[VX2] = v[VX3] = 0;
  v[BX1] = v[BX2] = v[BX3] = 0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  //int i, j, k;
  //int id1, id2;
  //double *x1 = grid->x[IDIR];
  //double *x2 = grid->x[JDIR];
  //double *x3 = grid->x[KDIR];
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid){
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0){
   B0[0] = 0.0;
   B0[1] = 1.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int call = 1;
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  int id1, id2;
  double Bext;
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  int call_x2_beg = 1;



  if (side == 0) {    /* -- check solution inside domain -- */
    TOT_LOOP(k,j,i){
    }
  }
  
  
  if (side == X2_END){  /* -- X2_END boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    //Bext = Bext0; //x2[j]/1.0e3 + Bext0;
    //if (g_time>=mt_start && g_time<=mt_end){
    //  if (x1[i]>mleft && x1[i]<mright){
    //    d->Vc[PRS][k][j][i] = PRS_1 +
    //        (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    //    d->Vc[RHO][k][j][i] = (T_f)*RHO_1;
    //    d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
    //    d->Vc[VX2][k][j][i] = -0.050*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
    //  }
    //}
    //else{
      d->Vc[BX2][k][j][i] = 0.0;
      d->Vc[RHO][k][j][i] = rho_min;
      d->Vc[PRS][k][j][i] = prs_min;
      d->Vc[VX2][k][j][i] = -0.8*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
      d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JEND-j-1][i];//1.0;
      d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JEND-j-1][i];//1.0;
    //d->flag[k][j][i]   |=      FLAG_INTERNAL_BOUNDARY;
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    //Bext = Bext0; //x2[j]/1.0e3 + Bext0;
    //if (g_time>=mt_start && g_time<=mt_end){
    //  if (x1[i]>mleft && x1[i]<mright){
    //    d->Vc[PRS][k][j][i] = PRS_0 +
    //        (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    //    d->Vc[RHO][k][j][i] = (T_f)*RHO_0;
    //    d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
    //  }
    //}
    //else{
      d->Vc[RHO][k][j][i] = rho_max;
      d->Vc[PRS][k][j][i] = prs_max;
      d->Vc[BX2][k][j][i] = 0.0;
      d->Vc[VX1][k][j][i] = d->Vc[VX1][k][2*JBEG-j+1][i];//1.0;
      d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
      d->Vc[BX1][k][j][i] = d->Vc[BX1][k][2*JBEG-j+1][i];//1.0;
    }
  }

}



#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  g[IDIR] = 0.0;
  //FILE *fptr;
  //fptr = fopen("GRAVITY0","a");
  //if (x1>0.9e4){
  //  fprintf(fptr,"%.5f %.9f\n",x2,interp_gra(x2));
  //  printf("%.5f %.9f\n",x2,interp_gra(x2));
  //}
  //fclose(fptr);
  g[JDIR] = 1.00*interp_gra(x2);
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
{
  return 0.0;
}
#endif




/* ********************************************************************* */
/* ********************************************************************* */

/* ADDED FUNCTIONS */


double interp_rho(double y){
  double rho_int;
  rho_int = 1.988984626e-07 + -1.073208823e-09*y +
  1.552630798e-11*pow(y,2) + 5.749153071e-13*pow(y,3) +
  7.226123210e-15*pow(y,4) + 4.994825132e-17*pow(y,5) +
  2.210451070e-19*pow(y,6) + 6.724332015e-22*pow(y,7) +
  1.462137306e-24*pow(y,8) + 2.321936803e-27*pow(y,9) +
  2.718388019e-30*pow(y,10) + 2.343529267e-33*pow(y,11) +
  1.469115589e-36*pow(y,12) + 6.511530032e-40*pow(y,13) +
  1.933542327e-43*pow(y,14) + 3.450190137e-47*pow(y,15) +
  2.796313533e-51*pow(y,16);  
  return rho_int/UNIT_DENSITY;
}


double interp_prs(double y){
  double prs_int;
  prs_int = 7.611740054e+04 + -5.459293398e+02*y +
  1.062505952e+00*pow(y,2) + -3.211659484e-02*pow(y,3) +
  -7.428620532e-04*pow(y,4) + -7.754892063e-06*pow(y,5) +
  -4.850736723e-08*pow(y,6) + -2.016973619e-10*pow(y,7) +
  -5.829629957e-13*pow(y,8) + -1.187332908e-15*pow(y,9) +
  -1.677769501e-18*pow(y,10) + -1.535992328e-21*pow(y,11) +
  -6.962603151e-25*pow(y,12) + 1.649177683e-28*pow(y,13) +
  3.712981400e-31*pow(y,14) + 6.750378344e-35*pow(y,15) +
  -1.407953614e-37*pow(y,16) + -6.172633205e-41*pow(y,17) +
  5.019330571e-44*pow(y,18) + 3.989473769e-47*pow(y,19) +
  -1.181344606e-50*pow(y,20) + -2.519249182e-53*pow(y,21) +
  -1.277193907e-56*pow(y,22) + -3.005252345e-60*pow(y,23) +
  -2.842338812e-64*pow(y,24);
  return prs_int/UNIT_PRESSURE;
}


double interp_gra(double y){
  double fun = -1.0*(274.0 - 0.0007873*y + 1.642e-09*pow(y,2));
  return fun*UNIT_DENSITY*UNIT_LENGTH*1e2/UNIT_PRESSURE;
}








