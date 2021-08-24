/* Attempt to evolve solar interior conditions from the PLUTO module. 
 * 
 * Solar function parameters are from Model S (Christensen-Dalsgaard) for 
 * the quiet Sun and from Cameron et al. 2011 and references therein for 
 * umbra and penumbra.
 *
 */

#include "pluto.h"

#define  x1_left      g_domBeg[0]
#define  x1_right     g_domEnd[0]
#define  x2_top       g_domBeg[1]
#define  x2_bottom    g_domEnd[1]
#define  mcenter      0.0e3      /* Center of magnetic field */
#define  mwidth       5.0e3     /* Width of magnetic field tube HWHM */
#define  mwidthp      10.0e3
#define  mleft        mcenter - 0.5*mwidth /* limits of the magnetic tube */
#define  mright       mcenter + 0.5*mwidth
#define  mleftp       mleft - 0.5*mwidthp
#define  mrightp      mright + 0.5*mwidthp
#define  mt_start     3600
#define  mt_end       3670


double UNIT_PRESSURE = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
double UNIT_MAGNETIC = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));

double interp_rho(double z);
double interp_prs(double z);
double interp_gra(double z);
double interp_rho_umbra(double z);
double interp_prs_umbra(double z);
double interp_rho_penumbra(double z);
double interp_prs_penumbra(double z);
double L_z(double z);  /* Logistic function */
double B_z(double z, double B0, double a, double a2);
double h_z(double z, double B0, double a, double a2, double h0);
double B_rz(double r, double z, double B0, double a, double a2, double h0);
double raised_cos(double z, double T, double beta, double z0);
double rho_min;  /* Values taken from the solar model at depth 0 and 10 Mm */
double rho_max;
double prs_min;
double prs_max;

//double mt_start = 600*60;  /* Start of magnetic field */
//double mt_end   = 600*60;

double B0 = 300;    /* KGauss Mag Field at center of tube */
double a  = 1250;   /* km */
//double a2 = 18400;  /* km */
double a2 = 1700;  /* km */
double h0 = mwidth;   /* HWHM of the gaussian at z=0 */
int call = 1;


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  int first_call = 1;

  /* Parameters for Paper simulation */
  double mag_inter;

  /* Call of upper and lower density, pressure variables */
  if (first_call){
    rho_min = interp_rho(x2_bottom);
    rho_max = interp_rho(x2_top);
    prs_min = interp_prs(x2_bottom);
    prs_max = interp_prs(x2_top);
    first_call = 0;
  }
  mag_inter = B_z(x2,B0,a,a2);
  //FILE *fptr;
  //fptr = fopen("raised.txt","a");
  //fprintf(fptr,"%.6f %.6f\n",x1,raised_cos(x1, 0.05e-3, 0.3, 0e3));
  //fprintf(fptr,"%.6f %.6f\n",x1,raised_cos(x1, 0.10e-3, 0.6, 15e3));
  //fprintf(fptr,"%.6f %.6f\n",x1,raised_cos(x1, 0.05e-3, 0.3, 30e3));
  //fprintf(fptr,"%.6f %.6f\n",x1,raised_cos(x1, 0.10e-3, 0.6, -15e3));
  //fprintf(fptr,"%.6f %.6f\n",x1,raised_cos(x1, 0.05e-3, 0.3, -30e3));
  //fclose(fptr);

  v[VX1] = v[VX2] = v[VX3] = 0.0;
  v[BX1] = v[BX3] = 0.0;
  //v[BX2] = 0.0;//B_rz(x1,x2,B0,a,a2,h0);


  if (g_time > mt_start){
    if (x1>mleft && x1<mright){
      v[RHO] = interp_rho_umbra(x2);
      v[PRS] = interp_prs_umbra(x2);
      v[BX2] = B_rz(x1,x2,B0,a,a2,h0);
    }
    else if ((x1>mleftp && x1<=mleft) || (x1>=mright && x1<mrightp)){
      v[RHO] = interp_rho_penumbra(x2);
      v[PRS] = interp_prs_penumbra(x2);
      v[BX2] = B_rz(x1,x2,B0,a,a2,h0);
    }
    else {
    //  v[BX2] = 0.0;
      v[RHO] = interp_rho(x2);
      v[PRS] = interp_prs(x2);
      //v[BX2] = B_z(x2,B0,a,a2);
      v[BX2] = B_rz(x1,x2,B0,a,a2,h0);
    }
  }
  else{
    v[RHO] = interp_rho(x2);
    v[PRS] = interp_prs(x2);
    //v[BX2] = B_z(x2,B0,a,a2);
    v[BX2] = B_rz(20e3,x2,B0,a,a2,h0);
  }


  //FILE *fptr;
  //fptr = fopen("PROFILE.txt","a");
  //fprintf(fptr,"%ld %.3f %.3f %.3f %.3f %.3f\n",g_stepNumber,x1/1e3,x2/1e3,
  //    v[RHO],v[PRS],mag_inter*UNIT_MAGNETIC);
  //fclose(fptr);

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid){
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
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  int id1, id2;
  double Bext;
  //static double ***rho_qs, ***prs_qs;
  //rho_qs = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  //prs_qs = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  int call_x2_beg = 1;

  if (side == 0) {    /* -- check solution inside domain -- */
    TOT_LOOP(k,j,i){
      if (g_time >= mt_start && g_time <= mt_end){
        //if (call){
        //  rho_qs[k][j][i] = d->Vc[RHO][k][j][i];
        //  prs_qs[k][j][i] = d->Vc[PRS][k][j][i];
        //  call = 0;
        //}
        if (x1[i]>mleft && x1[i]<mright){
          d->Vc[RHO][k][j][i] = interp_rho_umbra(x2[j]);//*raised_cos(x1[i], 0.05e-3, 0.3, 0e3);
          d->Vc[PRS][k][j][i] = interp_prs_umbra(x2[j]);//*raised_cos(x1[i], 0.05e-3, 0.3, 0e3);
          d->Vc[BX2][k][j][i] = B_rz(x1[i],x2[j],B0,a,a2,h0);
        }
        if (x1[i]>mleftp && x1[i]<=mleft){
          d->Vc[RHO][k][j][i] = interp_rho_penumbra(x2[j]);//*raised_cos(x1[i], 0.1e-3, 0.6, -15e3);
          d->Vc[PRS][k][j][i] = interp_prs_penumbra(x2[j]);//*raised_cos(x1[i], 0.1e-3, 0.6, -15e3);
          d->Vc[BX2][k][j][i] = B_rz(x1[i],x2[j],B0,a,a2,h0);
        }
        if (x1[i]>=mright && x1[i]<mrightp){
          d->Vc[RHO][k][j][i] = interp_rho_penumbra(x2[j]);//*raised_cos(x1[i], 0.1e-3, 0.6, 15e3);
          d->Vc[PRS][k][j][i] = interp_prs_penumbra(x2[j]);//*raised_cos(x1[i], 0.1e-3, 0.6, 15e3);
          d->Vc[BX2][k][j][i] = B_rz(x1[i],x2[j],B0,a,a2,h0);
        }
        //if (x1[i]>=mrightp){
        //  d->Vc[RHO][k][j][i] = rho_qs[k][j][i]*raised_cos(x1[i], 0.05e-3, 0.3, 30e3);
        //  d->Vc[PRS][k][j][i] = prs_qs[k][j][i]*raised_cos(x1[i], 0.05e-3, 0.3, 30e3);
        //}
        //if (x1[i]<=mleftp){
        //  d->Vc[RHO][k][j][i] = rho_qs[k][j][i]*raised_cos(x1[i], 0.05e-3, 0.3, -30e3);
        //  d->Vc[PRS][k][j][i] = prs_qs[k][j][i]*raised_cos(x1[i], 0.05e-3, 0.3, -30e3);
        //}
      }
    }
  }
  
  if (side == X2_END){  /* -- X2_END boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    //Bext = Bext0; //x2[j]/1.0e3 + Bext0;
    if (g_time>=mt_start){ // && g_time<=mt_end){
      if (x1[i]>mleft && x1[i]<mright){
        //  d->Vc[PRS][k][j][i] = PRS_1 +
        //      (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
        //  d->Vc[RHO][k][j][i] = (T_f)*RHO_1;
        //  d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
        //  d->Vc[VX2][k][j][i] = -0.050*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
        d->Vc[RHO][k][j][i] = interp_rho_umbra(x2[j]);
        d->Vc[PRS][k][j][i] = interp_prs_umbra(x2[j]);
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = B_rz(x1[i], x2[j], B0, a, a2, h0);
        //d->Vc[BX2][k][j][i] = B_rz(x1[i], x2[j], B0, a, a2, h0);
      }
      else if ((x1[i]>mleftp && x1[i]<=mleft) || (x1[i]>=mright && x1[i]<mrightp)){
        d->Vc[RHO][k][j][i] = interp_rho_penumbra(x2[j]);
        d->Vc[PRS][k][j][i] = interp_prs_penumbra(x2[j]);
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = B_rz(x1[i], x2[j], B0, a, a2, h0);
      }
      else{
        d->Vc[RHO][k][j][i] = rho_min;
        d->Vc[PRS][k][j][i] = prs_min;
        d->Vc[BX2][k][j][i] = B_rz(20e3, x2[j], B0, a, a2, h0);
        d->Vc[VX2][k][j][i] = -0.8*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
        d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JEND-j-1][i];//1.0;
        d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JEND-j-1][i];//1.0;
      }
    }
      else{
        d->Vc[RHO][k][j][i] = rho_min;
        d->Vc[PRS][k][j][i] = prs_min;
        d->Vc[BX2][k][j][i] = B_rz(20e3, x2[j], B0, a, a2, h0);
        d->Vc[VX2][k][j][i] = -0.8*d->Vc[VX2][k][2*JEND-j-1][i];//1.0;
        d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JEND-j-1][i];//1.0;
        d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JEND-j-1][i];//1.0;
      }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary: semi-permeable -- */
    BOX_LOOP(box,k,j,i){
    //Bext = Bext0; //x2[j]/1.0e3 + Bext0;
    if (g_time>=mt_start){ // && g_time<=mt_end){
      if (x1[i]>mleft && x1[i]<mright){
    //    d->Vc[PRS][k][j][i] = PRS_0 +
    //        (1.0/(2.0*pow(UNIT_PRESSURE,2)))*(pow(Bext,2)-pow(Bint,2));
    //    d->Vc[RHO][k][j][i] = (T_f)*RHO_0;
    //    d->Vc[BX2][k][j][i] = Bint/UNIT_PRESSURE;
        d->Vc[RHO][k][j][i] = interp_rho_umbra(x2[j]);
        d->Vc[PRS][k][j][i] = interp_prs_umbra(x2[j]);
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = B_rz(x1[i], x2[j], B0, a, a2, h0);
    //  d->Vc[BX2][k][j][i] = B_rz(x1[i], x2[j], B0, a, a2, h0);
      }
      else if ((x1[i]>mleftp && x1[i]<=mleft) || (x1[i]>=mright && x1[i]<mrightp)){
        d->Vc[RHO][k][j][i] = interp_rho_penumbra(x2[j]);
        d->Vc[PRS][k][j][i] = interp_prs_penumbra(x2[j]);
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[BX2][k][j][i] = B_rz(x1[i], x2[j], B0, a, a2, h0);
      }
      else{
        d->Vc[RHO][k][j][i] = rho_max;
        d->Vc[PRS][k][j][i] = prs_max;
        d->Vc[BX2][k][j][i] =  B_rz(20e3, x2[j], B0, a, a2, h0);
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][2*JBEG-j+1][i];//1.0;
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][2*JBEG-j+1][i];//1.0;
      }
    }
      else{
        d->Vc[RHO][k][j][i] = rho_max;
        d->Vc[PRS][k][j][i] = prs_max;
        d->Vc[BX2][k][j][i] = B_rz(20e3, x2[j], B0, a, a2, h0);
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][2*JBEG-j+1][i];//1.0;
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][2*JBEG-j+1][i];//1.0;
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][2*JBEG-j+1][i];//1.0;
      }
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


double interp_rho(double z){
  double rho_int;
  rho_int = 1.988984626e-07 + -1.073208823e-09*z +
  1.552630798e-11*pow(z,2) + 5.749153071e-13*pow(z,3) +
  7.226123210e-15*pow(z,4) + 4.994825132e-17*pow(z,5) +
  2.210451070e-19*pow(z,6) + 6.724332015e-22*pow(z,7) +
  1.462137306e-24*pow(z,8) + 2.321936803e-27*pow(z,9) +
  2.718388019e-30*pow(z,10) + 2.343529267e-33*pow(z,11) +
  1.469115589e-36*pow(z,12) + 6.511530032e-40*pow(z,13) +
  1.933542327e-43*pow(z,14) + 3.450190137e-47*pow(z,15) +
  2.796313533e-51*pow(z,16);  
  return rho_int/UNIT_DENSITY;
}


double interp_prs(double z){
  double prs_int;
  prs_int = 7.611740054e+04 + -5.459293398e+02*z +
  1.062505952e+00*pow(z,2) + -3.211659484e-02*pow(z,3) +
  -7.428620532e-04*pow(z,4) + -7.754892063e-06*pow(z,5) +
  -4.850736723e-08*pow(z,6) + -2.016973619e-10*pow(z,7) +
  -5.829629957e-13*pow(z,8) + -1.187332908e-15*pow(z,9) +
  -1.677769501e-18*pow(z,10) + -1.535992328e-21*pow(z,11) +
  -6.962603151e-25*pow(z,12) + 1.649177683e-28*pow(z,13) +
  3.712981400e-31*pow(z,14) + 6.750378344e-35*pow(z,15) +
  -1.407953614e-37*pow(z,16) + -6.172633205e-41*pow(z,17) +
  5.019330571e-44*pow(z,18) + 3.989473769e-47*pow(z,19) +
  -1.181344606e-50*pow(z,20) + -2.519249182e-53*pow(z,21) +
  -1.277193907e-56*pow(z,22) + -3.005252345e-60*pow(z,23) +
  -2.842338812e-64*pow(z,24);
  return prs_int/UNIT_PRESSURE;
}


double interp_gra(double z){
  double fun = -1.0*(274.0 - 0.0007873*z + 1.642e-09*pow(z,2));
  return fun*UNIT_DENSITY*UNIT_LENGTH*1e2/UNIT_PRESSURE;
}


double interp_prs_umbra(double z){
  double prs_int_umbra;
  if (z<-1500){
    prs_int_umbra = interp_prs(z);
  }
  else if (z>=-1500 && z<-650){
    prs_int_umbra = 9.928356382e+11 + 1.063497063e+10*z +
    4.867908794e+07*pow(z,2) + 1.215150159e+05*pow(z,3) +
    1.708623253e+02*pow(z,4) + 1.133238511e-01*pow(z,5) +
    -1.169989908e-05*pow(z,6) + -5.697791780e-08*pow(z,7) +
    -1.718321639e-12*pow(z,8) + 2.667088173e-14*pow(z,9) +
    -3.102709686e-18*pow(z,10) + -1.173677270e-20*pow(z,11) +
    5.109104006e-24*pow(z,12) + 3.643606371e-27*pow(z,13) +
    -3.902819219e-30*pow(z,14) + -1.823780990e-35*pow(z,15) +
    1.873737146e-36*pow(z,16) + -8.553251186e-40*pow(z,17) +
    -5.261642299e-43*pow(z,18) + 6.484094800e-46*pow(z,19) +
    1.846496755e-50*pow(z,20) + -3.418376779e-52*pow(z,21) +
    4.681352105e-56*pow(z,22) + 2.002453582e-58*pow(z,23) +
    8.624642263e-62*pow(z,24) + 1.180910004e-65*pow(z,25);
  }
  else if (z>=-650 && z<-300){
    prs_int_umbra = -1.529630685e+10 + -3.531461172e+08*z +
    -3.468431184e+06*pow(z,2) + -1.841168066e+04*pow(z,3) +
    -5.382545620e+01*pow(z,4) + -6.726404873e-02*pow(z,5) +
    5.071418647e-05*pow(z,6) + 2.054462681e-07*pow(z,7) +
    -7.503052218e-11*pow(z,8) + -5.163417471e-13*pow(z,9) +
    3.465088140e-16*pow(z,10) + 1.149772931e-18*pow(z,11) +
    -1.595670784e-21*pow(z,12) + -1.842041754e-24*pow(z,13) +
    5.545793051e-27*pow(z,14) + 6.015153389e-31*pow(z,15) +
    -1.525264286e-32*pow(z,16) + 8.055948383e-36*pow(z,17) +
    3.901695221e-38*pow(z,18) + -2.903539240e-41*pow(z,19) +
    -1.338646635e-43*pow(z,20) + -1.209814594e-46*pow(z,21) +
    -3.663098454e-50*pow(z,22);
  }
  else if (z>=-300 && z<=0){
    prs_int_umbra = 3.519885883e+03 + -1.248204492e+01*z +
    2.480949984e-01*pow(z,2) + 8.926816438e-03*pow(z,3) +
    1.565831576e-04*pow(z,4) + 1.437427086e-06*pow(z,5) +
    6.950269287e-09*pow(z,6) + 1.388634725e-11*pow(z,7) +
    -1.493947071e-14*pow(z,8) + -1.111677226e-16*pow(z,9) +
    -1.339625944e-19*pow(z,10);
  }
  else {
    prs_int_umbra = 3.515315743e+03 + -1.547491147e+01*z +
    1.073850747e-01*pow(z,2) + -5.431913317e-03*pow(z,3) +
    1.801440951e-04*pow(z,4) + -2.973644763e-06*pow(z,5) +
    2.689342954e-08*pow(z,6) + -1.319935307e-10*pow(z,7) +
    2.324248582e-13*pow(z,8) + 1.143490145e-15*pow(z,9) +
    -9.193043162e-18*pow(z,10) + 3.067081951e-20*pow(z,11) +
    -6.015829708e-23*pow(z,12) + 7.181110378e-26*pow(z,13) +
    -4.852040653e-29*pow(z,14) + 1.429180243e-32*pow(z,15);
  }
  return prs_int_umbra/UNIT_PRESSURE;
}
  

double interp_rho_umbra(double z){
  double rho_int_umbra;
  if (z<-1500){
    rho_int_umbra = interp_rho(z);
  }
  else if (z>=-1500 && z<-650){
    rho_int_umbra = -7.875116122e-02 + -8.646311569e-04*z +
    -4.152409908e-06*pow(z,2) + -1.127753159e-08*pow(z,3) +
    -1.845260089e-11*pow(z,4) + -1.716104371e-14*pow(z,5) +
    -5.634129757e-18*pow(z,6) + 4.810255312e-21*pow(z,7) +
    4.404715809e-24*pow(z,8) + -1.219043484e-27*pow(z,9) +
    -2.242566973e-30*pow(z,10) + 4.102288117e-34*pow(z,11) +
    1.095748058e-36*pow(z,12) + -1.551555740e-40*pow(z,13) +
    -5.573132571e-43*pow(z,14) + -8.396415025e-48*pow(z,15) +
    2.945438325e-49*pow(z,16) + 1.872063575e-52*pow(z,17) +
    4.977945273e-56*pow(z,18) + 5.124506244e-60*pow(z,19);
  }
  else if (z>=-650 && z<-300){
    rho_int_umbra = 3.973461619e-02 + 8.602732054e-04*z +
    7.904530113e-06*pow(z,2) + 3.908095950e-08*pow(z,3) +
    1.051536749e-10*pow(z,4) + 1.133421742e-13*pow(z,5) +
    -1.222588419e-16*pow(z,6) + -3.663807118e-19*pow(z,7) +
    1.901550999e-22*pow(z,8) + 9.229710148e-25*pow(z,9) +
    -6.205013985e-28*pow(z,10) + -2.165472134e-30*pow(z,11) +
    2.341691968e-33*pow(z,12) + 4.804982822e-36*pow(z,13) +
    -7.589801114e-39*pow(z,14) + -1.213014474e-41*pow(z,15) +
    1.930227313e-44*pow(z,16) + 5.165622239e-47*pow(z,17) +
    4.055698216e-50*pow(z,18) + 1.131111554e-53*pow(z,19);
  }
  else if (z>=-300 && z<=0){
    rho_int_umbra = 9.963924397e-09 + -3.452388999e-11*z +
    4.781406441e-13*pow(z,2) + 1.158748826e-14*pow(z,3) +
    3.189849738e-17*pow(z,4) + -1.991137404e-18*pow(z,5) +
    -3.014179800e-20*pow(z,6) + -2.036118022e-22*pow(z,7) +
    -7.434555541e-25*pow(z,8) + -1.423972505e-27*pow(z,9) +
    -1.121528914e-30*pow(z,10);
  }
  else {
    rho_int_umbra = 9.986203727e-09 + -4.828457151e-11*z +
    2.170428438e-13*pow(z,2) + 2.898805302e-15*pow(z,3) +
    -9.834498879e-17*pow(z,4) + 1.059846387e-18*pow(z,5) +
    -5.858421898e-21*pow(z,6) + 1.847207284e-23*pow(z,7) +
    -3.364678898e-26*pow(z,8) + 3.303882417e-29*pow(z,9) +
    -1.356382700e-32*pow(z,10);
  }
  return rho_int_umbra/UNIT_DENSITY;
}
  

double interp_rho_penumbra(double z){
  double rho_int_penumbra;
  if (z<-1500){
    rho_int_penumbra = interp_rho(z);
  }
  else if (z>=-1500 && z<-650){
    rho_int_penumbra = -3.195586756e-02 + -3.595423886e-04*z +
    -1.769214896e-06*pow(z,2) + -4.925086788e-09*pow(z,3) +
    -8.273064433e-12*pow(z,4) + -7.944723264e-15*pow(z,5) +
    -2.813475045e-18*pow(z,6) + 2.140741308e-21*pow(z,7) +
    2.140784000e-24*pow(z,8) + -5.067809966e-28*pow(z,9) +
    -1.093920844e-30*pow(z,10) + 1.622298633e-34*pow(z,11) +
    5.395867387e-37*pow(z,12) + -5.989450602e-41*pow(z,13) +
    -2.761456991e-43*pow(z,14) + -1.083992307e-47*pow(z,15) +
    1.454939401e-49*pow(z,16) + 9.497793393e-53*pow(z,17) +
    2.567304966e-56*pow(z,18) + 2.677558914e-60*pow(z,19);
  }
  else if (z>=-650 && z<-300){
    rho_int_penumbra = -1.202270769e-01 + -2.609620682e-03*z +
    -2.403664992e-05*pow(z,2) + -1.191359211e-07*pow(z,3) +
    -3.215015978e-10*pow(z,4) + -3.485852116e-13*pow(z,5) +
    3.712050895e-16*pow(z,6) + 1.125693614e-18*pow(z,7) +
    -5.723850862e-22*pow(z,8) + -2.840374097e-24*pow(z,9) +
    1.876857316e-27*pow(z,10) + 6.686296785e-30*pow(z,11) +
    -7.137623201e-33*pow(z,12) + -1.490067194e-35*pow(z,13) +
    2.325653630e-38*pow(z,14) + 3.771578201e-41*pow(z,15) +
    -5.927850930e-44*pow(z,16) + -1.601791250e-46*pow(z,17) +
    -1.263961965e-49*pow(z,18) + -3.541354287e-53*pow(z,19);
  }
  else if (z>=-300 && z<=0){
    rho_int_penumbra = 2.474287245e-08 + -3.540760494e-10*z +
    -8.870534822e-12*pow(z,2) + -4.190162026e-13*pow(z,3) +
    -8.798418653e-15*pow(z,4) + -1.053996732e-16*pow(z,5) +
    -7.556950340e-19*pow(z,6) + -3.308237961e-21*pow(z,7) +
    -8.652139357e-24*pow(z,8) + -1.240093362e-26*pow(z,9) +
    -7.482374484e-30*pow(z,10); 
  }
  else {
    rho_int_penumbra = 2.497745484e-08 + -2.549438712e-10*z +
    1.332017418e-12*pow(z,2) + -4.229985867e-15*pow(z,3) +
    6.217510346e-18*pow(z,4) + 1.767776887e-20*pow(z,5) +
    -1.545952141e-22*pow(z,6) + 5.269262215e-25*pow(z,7) +
    -1.007812919e-27*pow(z,8) + 1.041267391e-30*pow(z,9) +
    -4.508555499e-34*pow(z,10);
  }
  return rho_int_penumbra/UNIT_DENSITY;
}


double interp_prs_penumbra(double z){
  double prs_int_penumbra;
  if (z<-1500){
    prs_int_penumbra = interp_prs(z);
  }
  else if (z>=-1500 && z<-650){
    prs_int_penumbra = -1.186982683e+09 + -1.261913588e+07*z +
    -5.732444601e+04*pow(z,2) + -1.420424907e+02*pow(z,3) +
    -1.982377596e-01*pow(z,4) + -1.303010946e-04*pow(z,5) +
    1.407720401e-08*pow(z,6) + 6.541203621e-11*pow(z,7) +
    1.777780060e-15*pow(z,8) + -3.050231414e-17*pow(z,9) +
    3.551832406e-21*pow(z,10) + 1.340741917e-23*pow(z,11) +
    -5.793551028e-27*pow(z,12) + -4.197265893e-30*pow(z,13) +
    4.442035871e-33*pow(z,14) + 6.101690646e-38*pow(z,15) +
    -2.157041532e-39*pow(z,16) + 9.611494050e-43*pow(z,17) +
    6.227391291e-46*pow(z,18) + -7.445043235e-49*pow(z,19) +
    -3.258660437e-53*pow(z,20) + 3.989408720e-55*pow(z,21) +
    -4.940526457e-59*pow(z,22) + -2.352970053e-61*pow(z,23) +
    -1.032055767e-64*pow(z,24) + -1.434913257e-68*pow(z,25);
  }
  else if (z>=-650 && z<-300){
    prs_int_penumbra = -1.141387115e+11 + -2.585181258e+09*z +
    -2.487943215e+07*pow(z,2) + -1.291070825e+05*pow(z,3) +
    -3.668597523e+02*pow(z,4) + -4.340654976e-01*pow(z,5) +
    3.798964048e-04*pow(z,6) + 1.334956731e-06*pow(z,7) +
    -6.266692047e-10*pow(z,8) + -3.287406124e-12*pow(z,9) +
    2.584975581e-15*pow(z,10) + 7.001407258e-18*pow(z,11) +
    -1.082436615e-20*pow(z,12) + -1.020562391e-23*pow(z,13) +
    3.536805406e-26*pow(z,14) + -2.517961062e-31*pow(z,15) +
    -9.279684065e-32*pow(z,16) + 5.801390836e-35*pow(z,17) +
    2.303465491e-37*pow(z,18) + -1.930417254e-40*pow(z,19) +
    -7.946224818e-43*pow(z,20) + -6.940513914e-46*pow(z,21) +
    -2.047695695e-49*pow(z,22);
  }
  else if (z>=-300 && z<=0){
    prs_int_penumbra = 8.032410575e+03 + -1.200148880e+02*z +
   -2.949015319e+00*pow(z,2) + -1.421420409e-01*pow(z,3) +
   -2.967151519e-03*pow(z,4) + -3.555739776e-05*pow(z,5) +
   -2.547361298e-07*pow(z,6) + -1.112202447e-09*pow(z,7) +
   -2.893028276e-12*pow(z,8) + -4.110627772e-15*pow(z,9) +
   -2.449628677e-18*pow(z,10);
  }
  else {
    prs_int_penumbra = 8.110270073e+03 + -8.716346996e+01*z +
    5.729759744e-01*pow(z,2) + -5.812304255e-03*pow(z,3) +
    1.041559971e-04*pow(z,4) + -1.464363101e-06*pow(z,5) +
    1.356539704e-08*pow(z,6) + -8.442246204e-11*pow(z,7) +
    3.613949209e-13*pow(z,8) + -1.071883664e-15*pow(z,9) +
    2.172808876e-18*pow(z,10) + -2.872262347e-21*pow(z,11) +
    2.173905263e-24*pow(z,12) + -5.000265787e-28*pow(z,13) +
    -4.636909980e-31*pow(z,14) + 2.777092973e-34*pow(z,15);
  }
  return prs_int_penumbra/UNIT_PRESSURE;
}
    




/* Logistic function */
double L_z(double z){
  return 1/(1+exp(-1.0*z));
}

/* Magnetic field in the center of the disk */
double B_z(double z, double B0, double a, double a2){
  return B0*exp((-1.0*z/a)*L_z(-1.0*z/a)*L_z(z/a2))/UNIT_MAGNETIC;
}

/* Half With at Half Maximum (HWHM) of the gaussian function */
double h_z(double z, double B0, double a, double a2, double h0){
  return h0*sqrt(B0/(B_z(z,B0,a,a2)*UNIT_MAGNETIC));
}

/* Magnetic field as a function of depth and r */
double B_rz(double r, double z, double B0, double a, double a2, double h0){
  return B_z(z,B0,a,a2)*exp(-1.0*log(2)*r*r/pow(h_z(z,B0,a,a2,h0),2));
}



double raised_cos(double z, double T, double beta, double z0){
  double smooth;
  double b = beta;
  if (z>-z0 && z<=-(1.0+b)/(2.0*T)+z0){
    smooth = 0.0;}
  else if (z>-(1.0+b)/(2.0*T)+z0 && z<=-(1.0-b)/(2.0*T)+z0){
    smooth = (1.0/2.0)*(1.0+cos((CONST_PI*T/b)*(fabs(z-z0)-(1.0-b)/(2.0*T))));}
  else if (z>-(1.0-b)/(2.0*T)+z0 && z<=(1.0-b)/(2.0*T)+z0){
    smooth = 1.0;}
  else if (z>(1.0-b)/(2.0*T)+z0 && z<=(1.0+b)/(2.0*T)+z0){
    smooth = (1.0/2.0)*(1.0+cos((CONST_PI*T/b)*(fabs(z-z0)-(1.0-b)/(2.0*T))));}
  else{
    smooth = 0.0;}

  if (z0 == 30e3){
    if (z>=30e3){
    smooth = 1.0;}
  }
  else if (z0 == -30e3){
    if (z<=-30e3){
    smooth = 1.0;}
  }
  return smooth;
}

  




