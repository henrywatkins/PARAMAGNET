//#######################################################################
//
//                        GLOBAL VARIABLES
//
//      This file contains the globally defined variables for the
//      numerical program. It contains the parameters that are set
//      to define the situation/mesh/physics being solved.
//
//#######################################################################

#include "Global.h"
#include <cmath>

/* - - - - - - - - - - - - - - - - - - - - - - - -
 Define set of global variables common to all functions/routines
 - - - - - - - - - - - - - - - - - - - - - - - - */


 //###################################
 //         Physical state
 //###################################
 const double    n0=1e26;   //initial number density of electrons in m^-3
 const double    T0=20*1.6e-19;   //initial temperature of electrons in
 const int       Z=1;             //The ionization state of the plasma
 const double    period=0.5e-9;  //Time period of simulation in seconds
 const double    xDom=1.0e-3;     //Physical size of simulation in x in metres
 const double    yDom=1.0e-3;    //Physical size of simulation in y in metres
 const double    zDom=0.01e-3;     //Physical size of simulation in z in metres
 //###################################
 //      Collisional Parameters
 //###################################
 const double coulog=24-log(sqrt(n0/1e6)*1.6e-19/T0);
 const double tei0=5.3e39*pow(T0,1.5)/n0/Z/coulog;
 const double l0=tei0*sqrt(2*T0/9.1e-31);
 const double brag=1.33;
 const double mu=3.5e-14*l0*l0*n0;
 const double R=1758.2;
 //###################################
 //          Mesh variables
 //###################################
 const double    T=period/tei0; //Dimensionless time period of the simulation, units of the initial ei collision time
 const double    Lx=xDom/l0;    // x Length, in units of the initial collisional mean free path
 const double    Ly=yDom/l0;    // y length
 const double    Lz=zDom/l0;    // z length
 const int       m=3000;        //number of time steps
 const int       nx=128;         //number of space steps/mesh size in x
 const int       ny=128;      //mesh size in y
 const int       nz=4;         //mesh size in z
 const int      N=nx*ny*nz;       //Total number of mesh points
 const int      vars=17;       //the number of variables being solved
 const double   dt=T/(double)(m-1);
 const double   dx=Lx/(double)(nx-1);
 const double   dy=Ly/(double)(ny-1);
 const double   dz=Lz/(double)(nz-1);
 const double   alpha= dt/(2*dx);
 const double   beta= dt/(2*dy);
 const double   gammad= dt/(2*dz);
 const double   flim=0.06;
 //################################
 //      Laser variables
 //################################
 const double  wl=1e-6;             //wavelength
 const double  omega=2*3.14*3e8/wl;     //angular frequency
 const double  wn=2*3.14*l0/wl;     //normalised wavenumber
 const double  nc=omega*omega*3.14e-4/n0;    //normalised critical density
 const double  Intensity=1e18;             //max intensity in watts per square meter
 const double  Pc=0.0;//0.5*Intensity/T0/3e8/(omega*omega*3.14e-4);   //ponderomotive constant
 const double  IBc=Intensity/T0/3e8/(omega*omega*3.14e-4);               //IB constant
 const double  theta=1.0;//0.707;     ///cosine of angle between polarisation
 //################################
 //     Set Output variables
 //################################
const int tfrac=30;
const int nout=1;        //Denstiy
const int Vout=0;        //Velocity
const int Tout=1;        //Temperature
const int Eout=0;        //Electric Field
const int Bout=1;        //Magnetic Field
const int jout=0;        //Current
const int qout=1;        //Heat Flow
const int lout=1;
const int l2out=0;

//#################################
// Boundary Conditions
// Periodic=0
// Relfective=1
// fixed=2 - set to 0
// Outflow=3
//#################################
//x boundaries
const int  xbcon = 0;
//y boundaries
const int  ybcon = 0;
//z boundaries
const int  zbcon = 0;
