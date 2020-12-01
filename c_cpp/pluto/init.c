/* First attempt to evolve solar interior conditions from the PLUTO module. 
 * 
 *
 */


#include "pluto.h"

#define x_range   150  /* Pixels in x1-direction */
#define y_range   300  /* Pixels in x2-direction */

double gravity_vector[y_range];
double gravity_vector_in[y_range];
int    a          = 0;      /* Counter to initialize global gravity array */
double T_f        = 1.5;    /* Fraction of the temp. inside sunspot */
double Bext       = 0.1;    /* External magnetic field Gauss (constant)*/
double Bint       = 1000.0; /* Internal magnetic field Gauss (constant)*/
double radius     = 0.0e3;  /* Radius of the sunspot */
double press_unit = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));
int    radtopix(double x);
void   gravity_read_data(void);



/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  #if HAVE_ENERGY
  #endif

  v[TRC] = 0.0;
  /* Read gravity data from hydrostatic equilibrium. I set a counter just
   * to initialize gravity array */
  if (a==0){ 
    gravity_read_data(); 
  }

  //printf("gravity vector %.10f\n",gravity_vector[295]);
  //double rnd, randomN;
  //rnd = rand() % 200;
  //randomN = (rnd-100)/100;
  //v[VX1] = 40*randomN*exp(-(1e4-x2)/2e4*40)*exp(-(x1*x1)/(10e3*10e3));
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX3] = 0.0;
  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
  a +=1;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  int i, j, k;
  int id1, id2;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];
  
  id1 = InputDataOpen("density.dbl","grid.out"," ",0);
  id2 = InputDataOpen("pressure.dbl","grid.out"," ",0);
  TOT_LOOP(k,j,i){
    if (x1[i]>-radius && x1[i]<radius){ 
      d->Vc[BX2][k][j][i] = Bint/press_unit;
      d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = 1.00*(InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
          (1.0/(2.0*pow(press_unit,2)))*(pow(Bext,2)-pow(Bint,2)));
    }
    else {
      d->Vc[RHO][k][j][i] = InputDataInterpolate (id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]);
      d->Vc[BX2][k][j][i] = Bext/press_unit;
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
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  int id1, id2;
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  id1 = InputDataOpen("density.dbl","grid.out"," ",0);
  id2 = InputDataOpen("pressure.dbl","grid.out"," ",0);
  //printf("\ng_time %d\n",g_stepNumber);

  //TOT_LOOP(k,j,i){ /* To add a perturbation */
  //  if (g_stepNumber==10){
  //    double r2 = sqrt(pow(x1[i]+9e3,2)+pow(x2[j]-10e3,2));
  //    if (r2<5e1){
  //      //d->Vc[PRS][k][j][i] = 0.80;
  //      d->Vc[RHO][k][j][i] = 2.07e0;
  //      d->Vc[VX2][k][j][i] = -0.0;
  //      d->Vc[VX1][k][j][i] = 900.0;
  //    }
  //  }
  //}

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary: outflow -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary: outflow -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
  
  
  if (side == X2_END){  /* -- X2_END boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    if (x1[i]>-radius && x1[i]<radius){ 
      d->Vc[BX2][k][j][i] = Bint/press_unit;
      d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
          (1.0/(2.0*pow(press_unit,2)))*(pow(Bext,2)-pow(Bint,2));
      d->Vc[VX2][k][j][i] = -0.050*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
    }
    else{
      d->Vc[BX2][k][j][i] = d->Vc[BX2][k][2*JEND-j-1][i];//1.0;
      d->Vc[RHO][k][j][i] = InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate(id2,x1[i],x2[j],x3[k]);
      d->Vc[VX2][k][j][i] = -0.070*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
    }
      d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JEND-j-1][i];//1.0;
      d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JEND-j-1][i];//1.0;
      //d->flag[k][j][i]   |=      FLAG_INTERNAL_BOUNDARY;
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    if (x1[i]>-radius && x1[i]<radius){ 
      //d->Vc[BX2][k][j][i] = Bint/press_unit;
      //d->Vc[RHO][k][j][i] = (T_f)*InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      //d->Vc[PRS][k][j][i] = InputDataInterpolate (id2,x1[i],x2[j],x3[k]) +
          //(1.0/(2.0*pow(press_unit,2)))*(pow(Bext,2)-pow(Bint,2));
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i];
      d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i];
      d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
    }
    else{
      //d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i];
      //d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i];
      d->Vc[RHO][k][j][i] = InputDataInterpolate(id1,x1[i],x2[j],x3[k]);
      d->Vc[PRS][k][j][i] = InputDataInterpolate(id2,x1[i],x2[j],x3[k]);
      d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
    }
      d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JBEG-j+1][i];//1.0;
      d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JBEG-j+1][i];//1.0;
      d->Vc[BX2][k][j][i] =      d->Vc[BX2][k][2*JBEG-j+1][i];//1.0;
      //d->flag[k][j][i]   |=      FLAG_INTERNAL_BOUNDARY;
    }
  }
  InputDataClose(id1);
  InputDataClose(id2);

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}



#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
    int pix = radtopix(x2);
    if (x1>-radius && x1<radius){ 
      g[JDIR] = gravity_vector_in[pix];}
    else{
      g[JDIR] = gravity_vector[pix];}
    g[IDIR] = 0.0;
    g[KDIR] = 0.0;
    //printf("%f for pix %d\n",x2,pix);
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
{
  return 0.0;
}
#endif




/* ********************************************************************* */
/* ********************************************************************* */

/* ADDED FUNCTIONS BY THE CUTEST AND SMARTEST GUY */


/* This function reads the input of the solar model to get thermodynamical
 * variables for write the gravity into an array. The fraction 100.0/6.0 
 * accounts for the "zero" of the photosphere in the Solar Model data.
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
  id1 = InputDataOpen("density.dbl","grid.out"," ",0);
  id2 = InputDataOpen("pressure.dbl","grid.out"," ",0);
  for (int i=0; i<lim; i++){ /* 100/6 to correct for zero of radius */
    res            =       yran*(lim-i-1)/(double)lim + 100.0/6.0;
    res2           =       yran*i/(double)lim + 100.0/6.0;
    dens_new[i]    =       InputDataInterpolate(id1,0.0,res2,0.5);
    dens_new_in[i] = (T_f)*InputDataInterpolate(id1,0.0,res2,0.5);
    pres_new[i]    =       InputDataInterpolate(id2,0.0,res2,0.5);
    pres_new_in[i] =       InputDataInterpolate(id2,0.0,res2,0.5) +
        (1.0/(2.0*pow(press_unit,2)))*(pow(Bext,2)-pow(Bint,2));
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
    gravity_vector[i]    = -1.0/dens_new[i]*derivati[i];
    gravity_vector_in[i] = -1.0/dens_new_in[i]*derivati_in[i];
  }
  printf("grav %f  grav_in %f\n",gravity_vector[100],gravity_vector_in[100]);
  printf("End gravity array \n");
  //printf(  " dens %f  dens_in %f\n",dens_new[10],dens_new_in[10]);
  //printf(  " pres %f  pres_in %f\n",pres_new[10],pres_new_in[10]);
}
  

int    radtopix(double x){
  int     lim = 300;
  double yran = 1.0e4;
  double result = (double)lim*(x-100.0/6.0)/yran;
  return round(result);
}


