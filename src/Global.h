//####################################################
//
//              GLOBAL VARIABLES HEADER
//
//     This is just the header that is paired with
//            the Global.cpp data file
//
//####################################################

#ifndef Global_h
#define Global_h


//------------Reference-Values-----------------------------
//These are the initial values of the density and temperature that are used to make the variable dimensionless
extern const double n0,T0,flim,coulog0,l0,tei0,wl,wn,period,xDom,yDom,zDom,T,Lx,Ly,Lz,dt,dx,dy,dz,alpha,beta,gammad,brag,R,mu,omega,nc,Intensity,Pc,IBc, theta;
extern const int    Z,m,nx,ny,nz,N,vars;
//-----------Output-Definitions------------------------------
extern const int tfrac;

extern const int nout;
extern const int Vout;
extern const int Tout;
extern const int Eout;
extern const int Bout;
extern const int jout;
extern const int qout;
extern const int lout;
extern const int l2out;
//-----------Boundary-Conditions------------------------

//x boundaries
extern const int  xbcon;
//y boundaries
extern const int  ybcon;
//z boundaries
extern const int  zbcon;
//----------Tolerances-----------------------------------







#endif /* Global_h */
