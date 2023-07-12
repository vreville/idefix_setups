/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Polytropic spherical wind (compare w/ idefix) 
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
extern double ParkerWind(double);
extern double Get_PolyWind_NM(double, double, double);
extern double Get_Rcrit(double, double);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *********************************************************************** */
{
  double rstar,rhostar,pstar;
  double vesc,vkep,cs,PonRho;
  double rc, vc, vwind0, vwind;
  double OmegaZ; /* Rotation rate */
  double Bg[3]; /* Background field */
  double Sol[2]; /* Hydrostatic solution for photosphere */
  double gi; /* Polytropic index used for intialization */

  /**************************************/
  /***** Compute useful quantities *****/
  /**************************************/
  
  g_gamma = g_inputParam[GAMMA];
  rstar = g_inputParam[RSTAR];
  rhostar = g_inputParam[RHOSTAR];
  
  /* Set rotating frame */
  OmegaZ=g_inputParam[VROT_VESC]*sqrt(2.);
  #if ROTATING_FRAME
  g_OmegaZ = OmegaZ;
  OmegaZ=0.;
  #endif
  
  /* Sound speed. */
  cs = g_inputParam[CS_VESC]*sqrt(2.);
  gi=g_inputParam[POLYIDX];

  if(gi==1){
    /* Pressure to density ratio. */
    PonRho = cs*cs;
    /* Define critical radius = 0.25 vesc^2/cs^2 Rstar. */
    rc = 0.25 / (g_inputParam[CS_VESC]*g_inputParam[CS_VESC]);
    /* Calculate the wind speed at the surface of the star. */
    vwind0 = ParkerWind(1./rc) * cs;
    /* Get wind speed at current location. */
    vwind = ParkerWind(x1/rc) * cs;
  }
  else{
    /* Pressure to density ratio. */
    PonRho = cs*cs/gi;
    /* Critical radius */
    rc = 1.0/Get_Rcrit(gi,cs);    
    /* Define critical speed from r_crit */
    vc = sqrt(1.0*0.5/rc);
    /* Calculate the wind speed at the surface of the star. */
    vwind0 = Get_PolyWind_NM(gi,rstar,rc) * vc;
    /* Get wind speed at current location. */
    vwind = Get_PolyWind_NM(gi, x1 ,rc) * vc;
  }

  /****************************/
  /***** Set grid values. *****/
  /****************************/

  /* Set velocities. */
  us[VX1] = vwind;
  us[VX2] = 0.0;  
  us[VX3] = OmegaZ*x1*sin(x2);

  //printf("wind0 = %f, %f\n", vwind0, ParkerWind(1.0118401/rc)*cs);
  
  /* Set density based on mass continuity. */
  us[RHO] =rhostar * vwind0 /(vwind * x1 * x1);

  /* Set pressure */
  if (gi==1){
    us[PRS] = PonRho * us[RHO];
  }
  else{
    us[PRS] = PonRho*g_inputParam[RHOSTAR]*pow(us[RHO]/g_inputParam[RHOSTAR],gi);
  }
  
#if BACKGROUND_FIELD == YES
  us[BX1] = 0.0;
  us[BX2] = 0.0;
  us[BX3] = 0.0;
#else
  BackgroundField(x1, x2, x3, Bg);
  us[BX1] = Bg[0];
  us[BX2] = Bg[1];
  us[BX3] = Bg[2];
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
/* *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double *x1r, *x2r, *x3r;
  double *x1l, *x2l, *x3l;
  double var_int[NVAR];

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  x1r = grid->xr[IDIR]; x1l = grid->xl[IDIR];
  x2r = grid->xr[JDIR]; x2l = grid->xl[JDIR];
  x3r = grid->xr[KDIR]; x3l = grid->xl[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){}
  }

  if (side == X1_BEG){
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
	Init(var_int,x1[i],x2[j],x3[k]);
	NVAR_LOOP(nv){
	  d->Vc[nv][k][j][i]=var_int[nv];
	}
      }
    }
    else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
      BOX_LOOP(box,k,j,i){
	Init(var_int, x1[i], x2r[j], x3[k]);
	d->Vs[BX2s][k][j][i] = var_int[BX2];
      }
      #endif
    }
  }  
}

/* ************************************************************** */
void BackgroundField (double x1, double x2, double x3, double *B0)
/* 
 *
 * PURPOSE
 *
 *   Define the component of a static, curl-free background 
 *   magnetic field.
 *
 *
 * ARGUMENTS
 *
 *   x1, x2, x3  (IN)    coordinates
 *
 *   B0         (OUT)    vector component of the background field.
 *
 *
 **************************************************************** */
{  

  double mu = g_inputParam[VA_VESC]*sqrt(2.);
  double rss = g_inputParam[RSS];
  double rop,opFieldCoef;

  if(rss==0.0){
    B0[0] = 2.*mu*cos(x2)/(x1*x1*x1);
    B0[1] = mu*sin(x2)/(x1*x1*x1);
    B0[2] = 0.0;
  }
  else {
    if(x1 < rss){
      rop = x1;
      opFieldCoef = 1.0;
    }
    else {
      rop = rss;
      opFieldCoef=rss*rss/x1/x1;
    }

    B0[0] = 2*mu*(2*pow(rop,-3)+pow(rss,-3))/(pow(rss,-3)+2)*cos(x2)*opFieldCoef;
    B0[1] = 2*mu*(pow(rop,-3)-pow(rss,-3))/(pow(rss,-3)+2)*sin(x2)*opFieldCoef;
    B0[2] = 0.0;
  }  
}  

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/* *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/* *********************************************************************** */
{
  double phi;
  /*------------------------*/
  /*** star gravity field ***/
  /*------------------------*/

  /* Assumes Rsun = 1.0 */
  if(x1 < 1.0){

    phi = -pow(x1,2);
  } else{
    phi = -1.0/x1 ;
  }
  
  return phi;
}
#endif


/* **************************************************************** */
double Get_PolyWind_NM(double gamma, double r_sph, double r_crit)
/* 
 * Gives a solution for a polytropic 1D wind with Gamma between 1.004 and
 * 1.5. The output speed is normalized by the critical speed (speed at 
 * the critical point) and the input position x is as well normalized by
 * the critical radius : x=r/rc. Since rc is a function of Gamma, the 
 * user need the critical radius computed. This is done further 
 * in another function.
 * A Newton-Raphson method is used. 
 ****************************************************************** */
{

  double x,w,wref,f,fprime;

  if(r_sph < 0.5){
    return 0.001;
  }
  
  x=r_sph/r_crit;

  w=0.01;
  if (x >= 1){ w=10;}
  f=pow(w,gamma+1.)-pow(w,gamma-1.)*(4./x+(5.-3.*gamma)/(gamma-1))+2./(gamma-1)*pow(x,2.-2.*gamma); 

  while(fabs(f)>=1e-10){    
    fprime = (gamma+1.)*pow(w,gamma)-pow(w,gamma-2.)*(gamma-1.)*(4./x+(5.-3.*gamma)/(gamma-1.));
    wref=w;
    w=w-f/fprime;
    while(w <=0){w=(wref+w)/2;}
    f=pow(w,gamma+1.)-pow(w,gamma-1.)*(4./x+(5.-3.*gamma)/(gamma-1))+2./(gamma-1)*pow(x,2.-2.*gamma);
  }
  return w;
}
  
/* **************************************************************** */
double Get_Rcrit(double gamma, double v_sound)
/*
 *
 * Give the critical radius of a polytropic wind. 
 * Here rc=Rstar/x
 * Critical speed can be deduced as following : vc²=0.5*GMs/rc
 *
 * A Newton-Raphson method is used. 
 ****************************************************************** */
{
  
  double x,xref,f,fprime;
  double cs2;

  cs2=v_sound*v_sound;
  x=0.01;
  xref=10.0;
  f=pow(0.5/cs2*pow(x,3.-2.*gamma),(gamma+1.)/(gamma-1.))-(0.5/cs2*pow(x,3.-2.*gamma))*(4./x+(5.-3.*gamma)/(gamma-1.))+2./(gamma-1.)*pow(x,2-2.*gamma);
  
  while(fabs(f)>=1e-10){ // && abs(xref-x)>=1e-7){
    fprime=pow(0.5/cs2,(gamma+1.)/(gamma-1.))*(3.-2.*gamma)*(gamma+1.)/(gamma-1.)*pow(x,(4.-2.*gamma*gamma)/(gamma-1.))+(0.5/cs2*pow(x,3.-2.*gamma))*((4./x+(5.-3.*gamma)/(gamma-1.))*(2.*gamma-3.)*1./x+4./x/x)-4.*pow(x,1.-2.*gamma);
    xref=x;
    x=x-f/fprime;
    while(x<=0){x=(xref+x)/2.;}     
    f=pow(0.5/cs2*pow(x,3.-2.*gamma),(gamma+1.)/(gamma-1.))-(0.5/cs2*pow(x,3.-2.*gamma))*(4./x+(5.-3.*gamma)/(gamma-1.))+2./(gamma-1.)*pow(x,2-2.*gamma);
  }
  return x;
}

/**************************************************/
double ParkerWind(double x)
/*  Parker wind velocity in unit of iso sound speed
    x = radius / critical radius.
 **************************************************/
{
  double v, f;
  double vref;

  v = 1e-7;
  f = v*v-2*log(v)-4/x-4*log(x)+3;
  if (x>1) {v=10;}
  while (fabs(f) > 1e-10){
    vref = v;
    v = v - 0.5*f/(v-1/v);
    while (v < 0){
      v = (v + vref)/2;
    }
    f = v*v-2*log(v)-4/x-4*log(x)+3;
  }
  return v;
}

