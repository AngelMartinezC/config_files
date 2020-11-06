#include "pluto.h"

double gravity_vector[303];
int    radtopix(double x);
void   gravity_read_data(void);
int    a = 0;
double control;


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  //printf("\nHOLAAA %f \n",UNIT_VELOCITY);
  //printf("\nHOLAAA %f \n",sqrt(4*PI));
  double r;
  double B0 = 22.7/sqrt(4*PI);
  double B1 = 1.128941*B0;
  double P0 = B0*B0/(4*5/3);
  double P1 = B1*B1/(100*5/3);

  double r2;
  double x11 = x1-0.5;
  double x22 = x2-0.5;
  r2 = sqrt(x11*x11 + x22*x22);
  r = sqrt(x1*x1 + x2*x2);
  if (x1< 2.0e+3+4.0e+3 && x1>4.0e+3){ //6.0e+3
    v[BX2] = B0;
    v[PRS] = P0; //3*B0*B0/(4*5); //119.3662073189215;
    v[RHO] = 4.09961e-05;}
  else{
    v[BX2] = B1; //112.8941/sqrt(4*PI);
    v[PRS] = P1; //3*B1*B1/(100*5); //6.085325058980813;
    v[RHO] = 8.36e-6;}//1.70477956e-6;}
  
  double rnd, randomN;

  rnd = rand() % 200;//RandomNumber(0,100);
  randomN = (rnd-100)/100;
  //printf("%f\n",randomN);
 
  t0 = t0+1;
  //v[RHO] = 1.0e-4;
  //v[PRS] = 1.0;
  //v[BX2] = 0.0;
  //#if USE_RANDOM_PERTURBATION == YES
  //v[VX1] = 1*randomN*exp(-x2/2e4)*exp(-(x1-5e3)*(x1-5e3)/(10e3*10e3));
  v[VX1] = 0.3; //exp(-x2/2e4*35.0);
  //printf("%f\n",x2);
  printf("%f\n",*d);
  printf("\n");
  //v[VX1] = 0;
  //#else
    //v[VX1] = 0.0;
  //#endif
  //v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
  //v[PRS] = 1.0;
  #endif
  v[TRC] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD
  v[BX1] = 0.0;
  v[BX3] = 0.0;
  //v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
  #endif
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
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

  if (side == X1_END){  /* -- X1_END boundary -- */
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

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
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

  if (side == X2_END){  /* -- X2_END boundary -- */
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
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif


/* ********************************************************************* */
/* ********************************************************************* */
/* ********************************************************************* */

/* ADDED FUNCTIONS BY THE CUTEST AND SMARTEST GUY */

void   gravity_read_data(void){
  int    lim  = 300;
  double yran = 1.0e4;
  int    id1, id2;
  double pres_new[300];//[333][333][333];
  double dens_new[300];
  double radi_new[300];
  double derivati[300];
  double res;
  id1 = InputDataOpen("density.dbl","grid.out"," ",0);
  id2 = InputDataOpen("pressure.dbl","grid.out"," ",0);
  for (int i=0; i<lim; i++){
    res = yran*i/(double)lim + 100.0/6.0;
    pres_new[i] = InputDataInterpolate(id1,300.0,res,0.5);
    dens_new[i] = InputDataInterpolate(id2,300.0,res,0.5);
    radi_new[i] = res;
    printf("After \n");
  }
  InputDataClose(id1);
  InputDataClose(id2);
  double dp, dx;
  for (int i=0; i<lim; i++){
    if (i==0){
      dp = pres_new[i+1]-pres_new[i];
      dx = radi_new[i+1]-radi_new[i];
      derivati[i] = dp/dx;
    }
    else {
      dp = pres_new[i]-pres_new[i-1];
      dx = radi_new[i]-radi_new[i-1];
      derivati[i] = dp/dx;
    }
    gravity_vector[i] = -1.0/dens_new[i]*derivati[i];
  }
}
  

int    radtopix(double x){
  int     lim = 300;
  double yran = 1.0e4;
  double result = (double)lim*(x-100.0/6.0)/yran;
  return round(result);
}


