#include "idefix.hpp"
#include "setup.hpp"

real gammaGlob;
real cs_vescGlob;
real va_vescGlob;
real densityFloorGlob;
real ParkerWind(real x);

// User-defined boundaries
void UserDefBoundary(Hydro *hydro, int dir, BoundarySide side, real t) {

  DataBlock *data = hydro->data;
  
  if( (dir==IDIR) && (side == left)) {
    IdefixArray4D<real> Vc = hydro->Vc;
    IdefixArray4D<real> Vs = hydro->Vs;
    IdefixArray1D<real> x1 = data->x[IDIR];
    IdefixArray1D<real> x2 = data->x[JDIR];

    int ighost = data->nghost[IDIR];

    real cs=cs_vescGlob*sqrt(2.);
    real rc = 0.25 / (cs_vescGlob*cs_vescGlob);
    real vwind0 = ParkerWind(1./rc) * cs;
    real PonRho = cs*cs;
    real va_vesc = va_vescGlob;
    //real ParkerBC[ighost];

    hydro->boundary->BoundaryFor("UserDefBoundary", dir, side,
				 KOKKOS_LAMBDA (int k, int j, int i) {
				   real r = x1(i);
				   real th = x2(j);
				   real R=x1(i)*sin(x2(j));
				   real z=x1(i)*cos(x2(j));
				   //real vwind = ParkerBC[i] * cs;
				   real mu = va_vesc * sqrt(2.);
				  
				   Vc(RHO,k,j,i) = 1.0;//*vwind0/(vwind * r * r);
				   Vc(PRS,k,j,i) = PonRho * Vc(RHO, k, j, i);
				   Vc(VX1,k,j,i) = 0.0;//vwind;
				   Vc(VX2,k,j,i) = 0.0;
				   Vc(VX3,k,j,i) = 0.0;
				   Vc(BX1,k,j,i) = 2 * mu * cos(th)/(r*r*r);
				   Vc(BX2,k,j,i) = mu * sin(th)/(r*r*r);
				   Vc(BX3,k,j,i) = 0.0;
				   
				 });

    /*	idefix_for("UserDefBoundaryX1S",0,data->np_tot[KDIR],0,data->np_tot[JDIR]+1,0,ighost,
		   KOKKOS_LAMBDA (int k, int j, int i) {*/

    hydro->boundary->BoundaryForX2s("UserDefBoundaryBX2s", dir, side,
				    KOKKOS_LAMBDA(int k, int j, int i){

				      real r=data->x[IDIR](i);
				      real th=data->xl[JDIR](j);
				      real mu = va_vesc * sqrt(2.);

				      Vs(BX2s,k,j,i) = mu * sin(th)/(r*r*r);
		     
				    });
  }
}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {
  // Set the function for userdefboundary
  data.hydro->EnrollUserDefBoundary(&UserDefBoundary);
  //data.hydro->EnrollUserSourceTerm(&MySourceTerm);
  //data.hydro->EnrollInternalBoundary(&InternalBoundary);
  //output->EnrollUserDefVariables(&ComputeUserVars);
  
  gammaGlob=input.Get<real>("Hydro", "gamma", 0);
  cs_vescGlob=input.Get<real>("Setup", "cs_vesc", 0);
  va_vescGlob=input.Get<real>("Setup", "va_vesc", 0);

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
  // Create a host copy
  DataBlockHost d(data);

  // Make vector potential
  IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
  
  real x,y,z;
  
  real vphi,f,r,th;
  real cs=cs_vescGlob*sqrt(2.);
  real rc = 0.25 / (cs_vescGlob*cs_vescGlob);
  real vwind0 = ParkerWind(1./rc) * cs;

  printf("vwind0 = %f, %f\n", vwind0, ParkerWind(1.0118401/rc)*cs);
  real PonRho = cs*cs;
  real mu = va_vescGlob * sqrt(2.);    

  for(int k = 0; k < d.np_tot[KDIR] ; k++) {
    for(int j = 0; j < d.np_tot[JDIR] ; j++) {
      for(int i = 0; i < d.np_tot[IDIR] ; i++) {
	r=d.x[IDIR](i);
	th=d.x[JDIR](j);
	real R=r*sin(th);
	real z=r*cos(th);

	real vwind = ParkerWind(r/rc) * cs;

	d.Vc(RHO,k,j,i) = 1.0*vwind0/(vwind * r * r);
	d.Vc(PRS,k,j,i) = PonRho * d.Vc(RHO, k, j, i);
	d.Vc(VX1,k,j,i) = vwind;
	d.Vc(VX2,k,j,i) = 0.0;
	d.Vc(VX3,k,j,i) = 0.0;

	r=d.x[IDIR](i);
	th=d.xl[JDIR](j);
		  
	//d.Vs(BX1s, k , j, i) = 2 * mu * cos(th)/(r*r*r);
	//d.Vs(BX2s, k, j, i) = mu * sin(th)/(r*r*r);
	//d.Vs(BX3s, k, j, i) = 0.0;

	A(IDIR,k,j,i) = 0.0;
	A(JDIR,k,j,i) = 0.0;
	A(KDIR,k,j,i) = mu  * sin(th)/r/r;
      }
    }
  }

  // Make the field from the vector potential
  d.MakeVsFromAmag(A);
  
  // Send it all, if needed
  d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {
}

/**************************************************/
real ParkerWind(real x)
/*  Parker wind velocity in unit of iso sound speed                                                                                                                                                                                
    x = radius / critical radius.                                                                                                                                                                                                  
**************************************************/
{
  real v, f;
  real vref;

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
