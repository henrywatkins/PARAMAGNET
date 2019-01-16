//####################################################################
//
//                   TRANSPORT COEFFICIENTS HEADER
//
//   This header file is a pair of the TransportCoeff.cpp file,
//   it will declare these functions and allow the main program to
//   access their definitions at compile/linking time
//
//####################################################################

#ifndef TransportCoeff_h
#define TransportCoeff_h

extern const double alphavalues[8][10];

extern const double betavalues[8][11];

extern const double kappavalues[8][11];

double alpha_para(double chi, int Z);

double alpha_perp(double chi, int Z);

double alpha_wedg(double chi, int Z);

double beta_para(double chi, int Z);

double beta_perp(double chi, int Z);

double beta_wedg(double chi, int Z);

double kappa_para(double chi, int Z);

double kappa_perp(double chi, int Z);

double kappa_wedg(double chi, int Z);


double hall(int Z, double *vn, int position);


void setbdir(int position, double *vn, double *b);
void setcoeff(int Z, double chi,double *out,double *b,int q);


#endif /* TransportCoeff_h */
