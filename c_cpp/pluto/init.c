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


#define x_range   150  /* Pixels in x1-direction */
#define y_range   300  /* Pixels in x2-direction */
#define DEPTH   1.0e4  /* Depth from solar surface */

double gravity_vector[y_range];
double gravity_vector_in[y_range];
//double gravity_vector_new[300][150];//[NX2_TOT][NX3_TOT];
//int    counter  = 0;      /* Counter to initialize global gravity array */
double T_f      = 1.5;    /* Fraction of the temp. inside sunspot */
double Bext0    = 0.1;    /* External magnetic field Gauss (constant)*/
double Bint     = 1000.0; /* Internal magnetic field Gauss (constant)*/
double radius   = 0.0e3;  /* Radius of the perfectly circular umbra */
//double radius_step= 1.0e3;  /* Radius used for stepping. See SCHEME 5 */
//double radius1    = 0.0e3;  /* radius for different magnetic tubes */
//double radius2    = 0.0e3;
//double press_unit = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));
int  radtopix(double x);
void gravity_read_data(void);
//void gravity(const Data *d, double array[300][150]);



/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  /* Read just the gravity array. The control of the vector body force is
     made by reading the thermodynamical varaibles to calculate grad(P)/rho
     and (NOT YET) j\times B */
  static int first_call = 1;
  /* Read gravity data from hydrostatic equilibrium. */
  if (first_call){ 
    gravity_read_data(); 
  }
  v[VX1] = v[VX2] = v[VX3] = 0;
  v[BX1] = v[BX2] = 0;
  first_call = 0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  int i, j, k;
  int id1, id2;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  double Bext;
  
  id1 = InputDataOpen("density.dbl","grid.out","little",0,CENTER);
  id2 = InputDataOpen("pressure.dbl","grid.out","little",0,CENTER);
  TOT_LOOP(k,j,i){
    Bext = Bext0 ;//x2[j]/1.0e3 + Bext0;
    if (x1[i]>-radius && x1[i]<radius){ 
    //if (x1[i]>-radius2 && x1[i]<-radius1 || x1[i]>radius1 && x1[i]<radius2){ 
      d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
      d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
          (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    }
    else {
      d->Vc[RHO][k][j][i] = InputDataInterpolate (id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]);
      d->Vc[BX2][k][j][i] = Bext/UNIT_PRESSURE;
      //temp=(d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*KELVIN*0.672442+1.6e3;
    }
  }
  InputDataClose(id1);
  InputDataClose(id2);
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
  double RHO_0 = 721.5352397851;
  double PRS_0 = 600725.6887423340;
  double RHO_1 = 0.2737624535;
  double PRS_1 = 14.2569055828;
  //printf("grav bef %f\n",gravity_vector_new[50][0]);
  //gravity(d, gravity_vector_new);
  //printf("grav aft %f\n\n",gravity_vector_new[50][50]);
  //id1 = InputDataOpen("density.dbl","grid.out"," ",0,CENTER);
  //id2 = InputDataOpen("pressure.dbl","grid.out"," ",0,CENTER);
  //printf("\ng_time %d\n",g_stepNumber);

  if (side == 0) {    /* -- check solution inside domain -- */
    TOT_LOOP(k,j,i){
    Bext = Bext0; //x2[j]/1.0e3 + Bext0;
      //if (x1[i]>-radius && x1[i]<radius){ 
      //  //printf("%d  ",i);
      //  /* LIMITS ARE 69 and 83 */
      //  d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
      //  d->Vc[RHO][k][j][i] = T_f*d->Vc[RHO][k][j][85]; //(T_f);
      //  d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][85] + 
      //      (1.0/(2.0*pow(UNIT_PRESSURE,2)))*
      //      (pow(d->Vc[BX2][k][j][85],2)-
      //       pow(d->Vc[BX2][k][j][80],2));
      //  //d->Vc[PRS][k][j][i] = 1.085(InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
      //  //    (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2)));
      //}
      //else {
      //  d->Vc[RHO][k][j][i] = InputDataInterpolate (id1,x1[i],x2[j],x3[k]);
      //  d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]);
      //  d->Vc[BX2][k][j][i] = Bext/UNIT_PRESSURE;
      //  //temp=(d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*KELVIN*0.672442+1.6e3;
      //}
      //d->Vc[BX1][k][j][i] = 0.0;
      //counter += 1;
      //printf("counter %d\n",counter);
    }
  }
  
  
  if (side == X2_END){  /* -- X2_END boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    Bext = Bext0; //x2[j]/1.0e3 + Bext0;
    //if (g_time>=20880){
    //  if (x1[i]>-radius_step && x1[i]<radius_step){ 
    //    d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
    //    d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
    //    d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
    //        (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    //    d->Vc[VX2][k][j][i] = -0.050*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
    //  }
    //}
    if (x1[i]>-radius && x1[i]<radius){ 
    //if (x1[i]>-radius2 && x1[i]<-radius1 || x1[i]>radius1 && x1[i]<radius2){ 
      d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
      d->Vc[RHO][k][j][i] = (T_f)*RHO_1;
      d->Vc[PRS][k][j][i] = PRS_1 + 
          (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
      d->Vc[VX2][k][j][i] = -0.050*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
    }
    else{
      d->Vc[BX2][k][j][i] = Bext/UNIT_PRESSURE;
      //d->Vc[BX2][k][j][i] = d->Vc[BX2][k][2*JEND-j-1][i];//1.0;
      d->Vc[RHO][k][j][i] = RHO_1;
      d->Vc[PRS][k][j][i] = PRS_1;
      d->Vc[VX2][k][j][i] = -0.070*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
    }
      d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JEND-j-1][i];//1.0;
      d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JEND-j-1][i];//1.0;
      //d->flag[k][j][i]   |=      FLAG_INTERNAL_BOUNDARY;
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    Bext = Bext0; //x2[j]/1.0e3 + Bext0;
    //if (g_time>=20880){
    //  if (x1[i]>-radius_step && x1[i]<radius_step){ 
    //    d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
    //    d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
    //    d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
    //        (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    //  }
    //}
    if (x1[i]>-radius && x1[i]<radius){ 
    //if (x1[i]>-radius2 && x1[i]<-radius1 || x1[i]>radius1 && x1[i]<radius2){ 
      d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
      d->Vc[RHO][k][j][i] = (T_f)*RHO_0;
      d->Vc[PRS][k][j][i] = PRS_0 + 
          (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
      //d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i];
      //d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i];
      //d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
    }
    else{
      //d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i];
      //d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i];
      d->Vc[RHO][k][j][i] = RHO_0;
      d->Vc[PRS][k][j][i] = PRS_0;
      d->Vc[BX2][k][j][i] = Bext/UNIT_PRESSURE;
    }
      d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JBEG-j+1][i];//1.0;
      d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
      d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JBEG-j+1][i];//1.0;
      //d->flag[k][j][i]   |=      FLAG_INTERNAL_BOUNDARY;
    }
  }
  //InputDataClose(id1);
  //InputDataClose(id2);

}



#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  //int cass = 0;
  //gravity(v);
  int pix = radtopix(x2);
  double rho, rhoi;

  //printf("GRAV aft %f\n\n",gravity_vector_new[50][50]);
  //printf("x2  %f\n",x2);
  //rho = v[RHO];
  //if (call){
  ////  //printf("rho %.4f, step: %d,  x1: %11.4f,  x2: %10.4f,  JDIR: %d,  ss %d \n",rho,g_stepNumber,x1,x2,JDIR);
  //  printf("step: %d,  x1: %11.4f,  x2: %10.4f \n\n",g_stepNumber,x1,x2);
  ////  //counter += 1;
  //}
  //printf("%d %d\n",counter, g_stepNumber);

  //int   i, j, k, nv;
  //TOT_LOOP(k,j,i){
  //  //if (i==1){
  //  //  printf("j: %d\n",j);}
  //  //printf("i: %d,  j:%d,  k: %d,  rho: %f\n",i,j,k,rho);
  //  //printf("i: %d,  j:%d,  k: %d,  nv: %f\n",i,j,k,rho);
  //}
  ////if (x1>-radius2 && x1<-radius1 || x1>radius1 && x1<radius2){ 
  if (x1>-radius && x1<radius){ 
    g[JDIR] = gravity_vector_in[pix];}
  else{
    g[JDIR] = gravity_vector[pix];
  }
  g[IDIR] = 0.0;
  g[KDIR] = 0.0;
  //}
  //printf("%f for pix %d\n",x2,pix);
  //cass += 1;
  //counter += 1; 
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


/* This function reads the input of the solar model to get thermodynamical
 * variables for write the gravity into an array. The fraction 100.0/6.0 
 * accounts for the "zero" of the photosphere in the Solar Model data.
 * 
 * Check if really 100./6.0 is needed: Solved with rescaling of solar model
 * variables.
 *
 * Here I write the input double type files into an array in order to be 
 * constantly called instead of making the interpolation on every runtime.
 */
void   gravity_read_data(void){
  printf("Start gravity array\n");
  int    lim  = y_range;
  double yran = g_domEnd[JDIR];
  int    id1, id2;
  double pres_new[y_range], pres_new_in[y_range];
  double dens_new[y_range], dens_new_in[y_range];
  double radi_new[y_range], radi_new_in[y_range];
  double derivati[y_range], derivati_in[y_range];
  double res, res2;
  id1 = InputDataOpen("density.dbl","grid.out"," ",0,CENTER);
  id2 = InputDataOpen("pressure.dbl","grid.out"," ",0,CENTER);
  double Bext;
  for (int i=0; i<lim; i++){ /* 100/6 to correct for zero of radius */
    Bext = Bext0; //yran*i/(1.0e3*lim) + Bext0; /* DEPTH=?yran */
    res            =       yran*(lim-i-1)/(double)lim;// + 100.0/6.0;
    res2           =       yran*i/(double)lim;// + 100.0/6.0;
    dens_new[i]    =       InputDataInterpolate(id1,0.0,res2,0.5);
    dens_new_in[i] = (T_f)*InputDataInterpolate(id1,0.0,res2,0.5);
    pres_new[i]    =       InputDataInterpolate(id2,0.0,res2,0.5);
    pres_new_in[i] =       InputDataInterpolate(id2,0.0,res2,0.5) +
        (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    radi_new[i]    = res;
    radi_new_in[i] = res;
  }
  InputDataClose(id1);
  InputDataClose(id2);
  double dp, dp_in, dx, dx_in;
  for (int i=0; i<lim; i++){
    if (i==0){
      dp    = pres_new[i+1]    - pres_new[i];
      dp_in = pres_new_in[i+1] - pres_new_in[i];
      dx    = radi_new[i+1]    - radi_new[i];
      dx_in = radi_new_in[i+1] - radi_new_in[i];
      derivati[i]    = dp/dx;
      derivati_in[i] = dp_in/dx_in;
    }
    else {
      dp    = pres_new[i]    - pres_new[i-1];
      dp_in = pres_new_in[i] - pres_new_in[i-1];
      dx    = radi_new[i]    - radi_new[i-1];
      dx_in = radi_new_in[i] - radi_new_in[i-1];
      derivati[i]    = dp/dx;
      derivati_in[i] = dp_in/dx_in;
    }
    gravity_vector[i]    = -derivati[i]/dens_new[i];
    gravity_vector_in[i] = -derivati_in[i]/dens_new_in[i];
  }
  /* Check a random value */
  //printf("grav %f  grav_in %f\n",gravity_vector[100],gravity_vector_in[100]);
  printf("End gravity array \n");
  //printf(  " dens %f  dens_in %f\n",dens_new[10],dens_new_in[10]);
  //printf(  " pres %f  pres_in %f\n",pres_new[10],pres_new_in[10]);
}
  

int    radtopix(double x){
  int     lim = y_range;
  double yran = DEPTH;  /* 100/6 to correct for zero of radius photosphere */
  //double result = (double)lim*(x-100.0/6.0)/yran;
  double result = (double)lim*(x)/yran;
  return round(result);
}



void gravity(const Data *d, double array[300][150]){//, double grav[][NX2_TOT][NX3_TOT]){
  int i, j, k;
  double rho, prs;
  double dp;
  double dx = 3.33333333333;
  //rho = v[RHO];
  //printf("counter: %d  step: %d, NX2_TOT %d\n",counter,g_stepNumber,NX2_TOT);

  
  //double *dx, *dy;
  //dx = grid->dx[IDIR];
  //dy = grid->dx[JDIR];
  TOT_LOOP(k,j,i){
    rho = d->Vc[RHO][k][j][i];
  //  //prs = d->Vc[PRS][k][j][i];
  //  /* Derivative on the y-direction */
    if (j==0){
      dp = d->Vc[PRS][k][j+1][i] - d->Vc[PRS][k][j][i];}
    else {
      dp = d->Vc[PRS][k][j][i] - d->Vc[PRS][k][j-1][i];}
    array[j][i] = -dp/(dx*rho);
  //  //if (i==70){
  //  //  printf("dy %.5f, NX2_TOT %d,  j %d\n",dy[j],NX2_TOT,j);
  //  //}
  }
  //printf("%d %d %d\n",IBEG,JBEG,KBEG);
}




