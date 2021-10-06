/* Attempt to evolve solar interior conditions from the PLUTO module. 
 * 
 * Solar function parameters are from Model S (Christensen-Dalsgaard) for 
 * the quiet Sun and from Cameron et al. 2011 and references therein for 
 * umbra and penumbra.
 *
 */

#include "pluto.h"

#define  halt_frac_prs  0.1 //005
#define  halt_frac_rho  0.1
#define  mwidth         1.3e3     /* Width of magnetic field tube HWHM */


double UNIT_PRESSURE = UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY;
double UNIT_MAGNETIC;// = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));

double interp_rho(double z);
double interp_prs(double z);
double interp_gra(double z);
double halt_function(double x1, double x2, double r0, double width);


double h0 = mwidth;   /* HWHM of the gaussian at z=0 */

/* Limits to set constant prs, rho and B. This reduces wave reflections */
double z_top         = 500.0;
double z_bottom      = -8000; //-13000.0;


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  double z = x2;
  UNIT_MAGNETIC = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));

  v[VX1] = v[VX2] = 0.0;
  v[BX1] = v[BX2] = 0.0;

  if (x3>=z_top){
    v[VX3] = 0.0;
    v[RHO] = interp_rho(z_top)*(1-halt_frac_rho*
        halt_function(x1,z_top,0.0,mwidth));
    v[PRS] = interp_prs(z_top)*(1-halt_frac_prs*
        halt_function(x1,z_top,0.0,mwidth));
    v[BX3] = sqrt(2*halt_frac_prs*interp_prs(z_top)*
        halt_function(x1,z_top,0.0,mwidth)); 
  }
  else if (x3<=z_bottom){
    v[RHO] = interp_rho(z_bottom)*(1-halt_frac_rho*
        halt_function(x1,z_bottom,0.0,mwidth));
    v[PRS] = interp_prs(z_bottom)*(1-halt_frac_prs*
        halt_function(x1,z_bottom,0.0,mwidth));
    v[BX3] = sqrt(2*halt_frac_prs*interp_prs(z_bottom)*
        halt_function(x1,z_bottom,0.0,mwidth)); 
  }
  else{
    //v[BX2] = 1000*halt_function(x1,z,0.0,1.0e3)/UNIT_MAGNETIC;
    v[RHO] = interp_rho(x3)*(1-halt_frac_rho*
        halt_function(x1,z,0.0,mwidth));
    v[PRS] = interp_prs(x3)*(1-halt_frac_prs*
        halt_function(x1,z,0.0,mwidth));
    v[BX3] = sqrt(2*halt_frac_prs*interp_prs(x3)*
        halt_function(x1,z,0.0,mwidth)); 
  }
  //if (z<0) printf("%.4f  %.4f\n",z,v[RHO]);
  v[TRC]  = (x1<=1.3e3 ? 1.0:0.0);
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
void Analysis (const Data *d, Grid *grid){
  /* Energy analisis over the whole domain. The idea is to apply this
     setup into every averaged pixels */
  int i, j, k;
  double dV, vol, scrh;
  double Ekin, Eth_max, vx2, vy2, vz2;
  double *dx, *dy, *dz;

  /* ---- Set pointer shortcuts ---- */
  dx = grid->dx[IDIR];
  dy = grid->dx[JDIR];
  dz = grid->dx[KDIR];

  /* ---- Main loop ---- */
  Ekin = Eth_max = 0.0;

  DOM_LOOP(k,j,i){
    dV = dx[i]*dy[j]*dz[k];

    /* Cell volume (Cartesian coordinates) */
    vx2 = d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i]; /* x-velocity squared */
    vy2 = d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i]; /* y-velocity squared */
    vz2 = d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]; /* z-velocity squared */

    scrh  = 0.5*d->Vc[RHO][k][j][i]*(vx2 + vy2 + vz2); /* cell kinetic energy */
    Ekin += scrh*dV;
    scrh  = d->Vc[PRS][k][j][i]/(g_gamma - 1.0); /* cell internal energy */
    Eth_max = MAX(Eth_max, scrh);
  }
  vol  = g_domEnd[IDIR] - g_domBeg[IDIR]; /* Compute total domain volume */
  vol *= g_domEnd[JDIR] - g_domBeg[JDIR];
  vol *= g_domEnd[KDIR] - g_domBeg[KDIR];

  Ekin /= vol; /* Compute kinetic energy average */
 
  /* ---- Parallel data reduction ---- */
  #ifdef PARALLEL
    MPI_Allreduce (&Ekin, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Ekin = scrh;
    MPI_Allreduce (&Eth_max, &scrh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    Eth_max = scrh;
    MPI_Barrier (MPI_COMM_WORLD);
  #endif

  /* ---- Write ascii file "averages.dat" to disk ---- */
  if (prank == 0){
    char fname[512];
    static double tpos = -1.0;
    FILE *fp;
    sprintf (fname, "%s/averages.dat",RuntimeGet()->output_dir);
    if (g_stepNumber == 0){ /* Open for writing only when we’re starting */
      fp = fopen(fname,"w"); /*
      from beginning */
      fprintf (fp,"# %7s %12s %12s %12s\n", "t", "dt", "<Ekin>","Max(Eth)");
    } else{
      /* Append if this is not step 0 */
      if (tpos < 0.0){ /* Obtain time coordinate of to last written row */
      char
      sline[512];
      fp = fopen(fname,"r");
      while (fgets(sline, 512, fp)) {}
      sscanf(sline, "%lf\n",&tpos); /* tpos = time of the last written row */
      fclose(fp);
    }
    fp = fopen(fname,"a");
  }
  if (g_time > tpos){ /* Write if current time if > tpos */
    fprintf (fp, "%12.6e %12.6e %12.6e %12.6e \n",g_time, g_dt,Ekin, Eth_max);
  }
  fclose(fp);
  }
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
  UNIT_MAGNETIC = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double r, r2;
  double pikes, sigma;
  double rnum1, rnum2, rnum3;
  double T, nu;
  x1 = grid->x[IDIR];  //[i];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];  //[k];
  double z1, z0;

  if (side == 0) {    /* -- check solution inside domain -- */

    DOM_LOOP(k,j,i){
      if (x3[k]>=z_top){
        d->Vc[RHO][k][j][i] = interp_rho(z_top)*(1-halt_frac_rho*
            halt_function(x1[i],z_top,0.0,mwidth));
        d->Vc[PRS][k][j][i] = interp_prs(z_top)*(1-halt_frac_prs*
            halt_function(x1[i],z_top,0.0,mwidth));
        d->Vc[BX3][k][j][i] = sqrt(2*halt_frac_prs*interp_prs(z_top)*
            halt_function(x1[i],z_top,0.0,1.0e3)); 
        d->Vc[VX3][k][j][i] = 0.0;
      }

      else if (x3[k]<=z_bottom){
        d->Vc[RHO][k][j][i] = interp_rho(z_bottom)*(1-halt_frac_rho*
            halt_function(x1[i],z_bottom,0.0,mwidth));
        d->Vc[PRS][k][j][i] = interp_prs(z_bottom)*(1-halt_frac_prs*
            halt_function(x1[i],z_bottom,0.0,mwidth));
        d->Vc[BX3][k][j][i] = sqrt(2*halt_frac_prs*interp_prs(z_bottom)*
            halt_function(x1[i],z_bottom,0.0,mwidth)); 
        d->Vc[VX3][k][j][i] = 0.0;
      }
    }

    if (g_time>=0.0*10){// && g_time<=20*60){
      sigma = 0.005;
      T = 6670;
      nu = 1.0/(2*T);
      rnum1 = 0.7e3; //(rand() % (int)16e3) - 8e3;
      rnum2 = -1500; //(rand() % (int)0.5e3) - 5.3e3; /* Depth */
      //rnum3 = 0.0005; //1.0 + ((rand() % (int)1) - 1.0)/10.0; /* Depth */
      rnum3 = 0.0001; //1.0 + ((rand() % (int)1) - 1.0)/10.0; /* Depth */
      pikes = exp(-pow(5*sin(2*CONST_PI*nu*(g_time)),2)/pow(sigma,2))*rnum3;
      //printf("pikes %.4f\n",pikes);
      DOM_LOOP(k,j,i){
        //r  = sqrt(pow(x1[i]-rnum1,2)   + 0.05*pow(x3[k]-rnum2,2));
        r  = sqrt(pow(x1[i]*cos(x2[j])-rnum1,2) + pow(x1[i]*sin(x2[j]),2) + 0.05*pow(x3[k]-rnum2,2));
        if (r<195){
          //d->Vc[PRS][k][j][i] += interp_prs(rnum2)*pikes; 
          d->Vc[BX3][k][j][i] += sqrt(2*halt_frac_prs*interp_prs(rnum2)*
              halt_function(x1[i],rnum2,0.0,mwidth))*pikes;
          d->Vc[VX3][k][j][i] += -0.00001;
        }
      }
    }
  }
  

  //if (side == X3_BEG){  /* -- X2_BEG boundary: semi-permeable -- */
  //  X3_BEG_LOOP(k,j,i) {  
  //    d->Vc[PRS][k][j][i] =      d->Vc[PRS][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[RHO][k][j][i] =      d->Vc[RHO][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[VX1][k][j][i] =      d->Vc[VX1][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[VX2][k][j][i] = -1.0*d->Vc[VX2][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[VX3][k][j][i] =      d->Vc[VX3][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[BX1][k][j][i] =      d->Vc[BX1][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[BX2][k][j][i] =      d->Vc[BX2][k][2*JBEG-j-1][i];//1.0;
  //    d->Vc[BX3][k][j][i] =      d->Vc[BX3][k][2*JBEG-j-1][i];//1.0;
  //  }
  //}
  //else if (side == X2_END){  /* -- X2_BEG boundary: semi-permeable -- */
  //  X2_END_LOOP(k,j,i) {  
  //  }
  //}

}



#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  g[IDIR] = 0.0;
  //fptr = fopen("GRAVITY0","a");
  //if (x1>0.9e4){
  //  fprintf(fptr,"%.5f %.9f\n",x2,interp_gra(x2));
  //  printf("%.5f %.9f\n",x2,interp_gra(x2));
  //}
  //fclose(fptr);
  g[JDIR] = 0.0;
  g[KDIR] = 1.00*interp_gra(x3);
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
  if (z>=-1500 && z<=0){
    rho_int = 1.988984626e-07 + -1.073208823e-09*z +
    1.552630798e-11*pow(z,2) + 5.749153071e-13*pow(z,3) +
    7.226123210e-15*pow(z,4) + 4.994825132e-17*pow(z,5) +
    2.210451070e-19*pow(z,6) + 6.724332015e-22*pow(z,7) +
    1.462137306e-24*pow(z,8) + 2.321936803e-27*pow(z,9) +
    2.718388019e-30*pow(z,10) + 2.343529267e-33*pow(z,11) +
    1.469115589e-36*pow(z,12) + 6.511530032e-40*pow(z,13) +
    1.933542327e-43*pow(z,14) + 3.450190137e-47*pow(z,15) +
    2.796313533e-51*pow(z,16);  
  }
  else if (z>0){
    rho_int = 1.998575400e-07 + -1.141330474e-09*z +
    1.405553852e-12*pow(z,2) + 2.641820862e-14*pow(z,3) +
    -1.125662086e-15*pow(z,4) + 2.466940533e-17*pow(z,5) +
    -3.097585201e-19*pow(z,6) + 2.503689907e-21*pow(z,7) +
    -1.392260630e-23*pow(z,8) + 5.507524992e-26*pow(z,9) +
    -1.568564535e-28*pow(z,10) + 3.199879633e-31*pow(z,11) +
    -4.566513080e-34*pow(z,12) + 4.331552831e-37*pow(z,13) +
    -2.453883829e-40*pow(z,14) + 6.283097143e-44*pow(z,15);
  }
  else if (z>-10000 && z<-1500){
    rho_int = -1.429784133e-05 + -4.021439051e-08*z +
    -3.767142430e-11*pow(z,2) + -7.077268240e-15*pow(z,3) +
    1.907801314e-17*pow(z,4) + 2.046690339e-20*pow(z,5) +
    1.092018763e-23*pow(z,6) + 3.638924192e-27*pow(z,7) +
    7.917976094e-31*pow(z,8) + 1.078861446e-34*pow(z,9) +
    6.939700989e-39*pow(z,10) + -4.159171351e-43*pow(z,11) +
    -1.438512650e-46*pow(z,12) + -1.549341076e-50*pow(z,13) +
    -9.159478994e-55*pow(z,14) + -2.970726413e-59*pow(z,15) +
    -4.155657549e-64*pow(z,16);
  }
  else if (z<=-10000){
    rho_int = -1.821979478e-01 + -1.158232305e-04*z +
    -3.274046534e-08*pow(z,2) + -5.420969401e-12*pow(z,3) +
    -5.817040917e-16*pow(z,4) + -4.227588907e-20*pow(z,5) +
    -2.107368013e-24*pow(z,6) + -7.115246845e-29*pow(z,7) +
    -1.557591565e-33*pow(z,8) + -1.996840308e-38*pow(z,9) +
    -1.138859985e-43*pow(z,10);
  }
  return rho_int/UNIT_DENSITY;
}


double interp_prs(double z){
  double prs_int;
  if (z>=-1500 && z<=0){
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
  }
  else if (z>0){
    prs_int = 7.609608761e+04 + -5.464792423e+02*z +
    1.374769105e+00*pow(z,2) + 7.920870429e-03*pow(z,3) +
    -2.564518093e-04*pow(z,4) + 4.454516117e-06*pow(z,5) +
    -5.244194845e-08*pow(z,6) + 4.266910248e-10*pow(z,7) +
    -2.452404978e-12*pow(z,8) + 1.011067030e-14*pow(z,9) +
    -3.004004898e-17*pow(z,10) + 6.380668296e-20*pow(z,11) +
    -9.453039879e-23*pow(z,12) + 9.279388478e-26*pow(z,13) +
    -5.423937333e-29*pow(z,14) + 1.428992380e-32*pow(z,15);
  }
  else if (z>-10000 && z<-1500){
    prs_int = 6.520571918e+07 + 2.447670723e+05*z +
    4.139831924e+02*pow(z,2) + 4.108812326e-01*pow(z,3) +
    2.702772915e-04*pow(z,4) + 1.226536873e-07*pow(z,5) +
    3.956676322e-11*pow(z,6) + 9.047798391e-15*pow(z,7) +
    1.424963548e-18*pow(z,8) + 1.391275422e-22*pow(z,9) +
    4.430506720e-27*pow(z,10) + -9.057682934e-31*pow(z,11) +
    -1.625523481e-34*pow(z,12) + -1.344203189e-38*pow(z,13) +
    -6.461427454e-43*pow(z,14) + -1.737974124e-47*pow(z,15) +
    -2.032311785e-52*pow(z,16);
  }
  else if (z<=-10000){
    prs_int = -9.225197855e+10 + -5.305395316e+07*z +
    -1.335816699e+04*pow(z,2) + -1.922869720e+00*pow(z,3) +
    -1.720445009e-04*pow(z,4) + -9.957023566e-09*pow(z,5) +
    -3.650899521e-13*pow(z,6) + -7.799710669e-18*pow(z,7) +
    -7.071511615e-23*pow(z,8) + 3.747156554e-28*pow(z,9) +
    9.379712599e-33*pow(z,10);
  }
  return prs_int/UNIT_PRESSURE;
}


double interp_gra(double z){
  double fun = -1.0*(274.0 - 0.0007873*z + 1.642e-09*pow(z,2));
  return fun*UNIT_DENSITY*UNIT_LENGTH*1e2/UNIT_PRESSURE;
}


double halt_function(double x1, double x2, double r0, double width){
  double func1, func2;
  double b = 20;
  double r = sqrt((x1-r0)*(x1-r0)); // + (x2-r0)*(x2-r0));
  /* before with 0.01 */
  //func1 = 1.0/(1.0+exp(0.09*(r-width)));
  func1 = 0.5*(1.0-tanh((r/width-1)*b));
  //func1 = 1.0/(1.0+exp(0.006*(r-width)))*1.0/(1.0+exp(0.003*(-x2-5e3)));
  return func1;
}


/* This is previous work in order to fit pressure and density from the
   magnetic field. This field has an inclination of ~pi/4 at z=0, and is
   modeled according to logistic functions */

///* Logistic function  */
//double L_z(double z){
//  return 1/(1+exp(-1.0*z));
//}
//double L_z3(double z){
//  double k = 2.0;
//  return 1/(1+exp(-1.0*k*(z+k)));
//}
//
///* Magnetic field in the center of the disk */
//double B_z(double z, double B0, double a, double a2){
//  UNIT_MAGNETIC = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));
//  return B0*L_z3(-1.0*z/a0)*exp((-1.0*z/a)*L_z(-1.0*z/a)*L_z(z/a2))/UNIT_MAGNETIC;
//}
//
///* Half With at Half Maximum (HWHM) of the gaussian function */
//double h_z(double z, double B0, double a, double a2, double h0){
//  UNIT_MAGNETIC = sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2));
//  return h0*sqrt(B0/(B_z(z,B0,a,a2)*UNIT_MAGNETIC));
//}
//
///* Magnetic field as a function of depth and r */
//double B_rz(double r, double z, double B0, double a, double a2, double h0){
//  return B_z(z,B0,a,a2)*exp(-1.0*log(2)*r*r/pow(h_z(z,B0,a,a2,h0),2));
//}
//
//double interp_rho_new(double r, double z){
//  double dB_dz = B0*(exp(-z/(a*(exp(z/a) + 1)*(exp(-z/a2) + 1)))*(-1*(exp(z/a2)*(exp(z/a0) + 1)*(a*(exp(z/a) + 1)*(a2*exp(z/a2) + a2 + z) - a2*z*exp(z/a)*(exp(z/a2) + 1)))/(a*a*a2*pow((exp(z/a)+1),2)*pow((exp(z/a2) + 1),2)) - exp(z/a0)/a0))/pow((exp(z/a0) + 1),2);
//
//  double dh_dz = h0*((exp(z/a0) + 1)*exp(z/(a*(exp(z/a) + 1)*(exp(-z/a2) + 1)))*(-(z*exp(z/a))/(a*a*pow((exp(z/a) + 1),2)*(exp(-z/a2) + 1)) + (z*exp(-z/a2))/(a*a2* (exp(z/a) + 1)*pow((exp(-z/a2) + 1),2)) + 1/(a*(exp(z/a) + 1)*(exp(-z/a2) + 1))) + exp(z/(a*(exp(z/a) + 1)*(exp(-z/a2) + 1)) + z/a0)/a0)/(2* sqrt((exp(z/a0) + 1)*exp(z/(a*(exp(z/a) + 1)*(exp(-z/a2) + 1)))));
//
//  double expon = exp(-1*2*log(2)*r*r/pow(h_z(z,B0,a,a2,h0),2));
//
//  double rho_res0 = interp_rho(z);
//
//  double rho_res1 = interp_rho(z) - 1e-2*(1.0/(2*interp_gra(z))*expon*B_z(z,B0,a,a2)*( 2*dB_dz + B_z(z,B0,a,a2)*(4*log(2)*r*r/pow(h_z(z,B0,a,a2,h0),3))*dh_dz ));
//
//  return rho_res1;
//}


