/* PLUTO code to evaluate the functionality of openmpi v.2.1.6
 * 
 * Here I evolve a drop with greater density than sourronging environment.
 */

#include "pluto.h"


int frame=0;
/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  v[RHO] = 1.0;
  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;
  #if HAVE_ENERGY
  v[PRS] = 1.0;
  #endif
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
  //printf("frame %d\n",frame);
  frame += 1;
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  printf("Step time: %d.  Code time: %f\n",g_stepNumber,g_time);

  if (g_time > 50 && g_time < 60){
    TOT_LOOP(k,j,i){ /* To add a perturbation */
    double r2 = sqrt(pow(x1[i]-800,2)+pow(x2[j]-550,2));
    if (r2<10e1){
      d->Vc[PRS][k][j][i] = 1.25;
      d->Vc[RHO][k][j][i] = 1.35;
      d->Vc[VX1][k][j][i] = -3.5;
      d->Vc[VX2][k][j][i] = -1.5;
      }
    }
  }
  if (g_time > 50 && g_time < 60){
    TOT_LOOP(k,j,i){ /* To add a perturbation */
    double r2 = sqrt(pow(x1[i]-500,2)+pow(x2[j]-500,2));
    if (r2<10e1){
      d->Vc[PRS][k][j][i] = 1.25;
      d->Vc[RHO][k][j][i] = 1.35;
      d->Vc[VX1][k][j][i] = -3.0;
      d->Vc[VX2][k][j][i] = -1.0;
      }
    }
  }
  else if (g_time > 80 && g_time < 110){
    TOT_LOOP(k,j,i){ /* To add a perturbation */
    double r3 = sqrt(pow(x1[i]- 70,2)+pow(x2[j]-200,2));
    if (r3<10e1){
      d->Vc[PRS][k][j][i] = 1.25;
      d->Vc[RHO][k][j][i] = 1.55;
      d->Vc[VX1][k][j][i] = 4.0;
      d->Vc[VX2][k][j][i] = 3.0;
      }
    }
  }



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
