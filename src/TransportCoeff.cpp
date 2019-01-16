//######################################################################################################
//
//                                      TRANSPORT COEFFICIENTS
//
//   These functions evaluate the dimensionless Transport Coefficients for the case of the magnetized
//   plasma in terms of the dimensionless hall parameter chi. The dependence on chi and values of the
//   coefficients were derived by Epperlein and Haines (Phys. Fluids, 1986). The coefficient values
//                                 are functions of the ionization.
//
//######################################################################################################

#include <math.h>
#include "TransportCoeff.h"
#include "Global.h"

using namespace std;

//Since the coefficients are functions of Z their values for different Z will be set into a database
//and pulled out when called in the functions below. The values were only calculated for Z of 1-8,10,
//12,14,20,30,60. The ones given here will be for 1-8. The row number will be Z, the column the coefficient type


const double alphavalues[8][10] = {{0.5061,1.37,3.03,2.77,6.72,2.66e2,2.53,3.28e3,3.46e3,3.66e2},{0.4295,1.58,3.21,2.78,6.70,4.91e2,2.53,3.72e3,6.69e3,5.54e2},{0.3950,1.68,3.17,2.78,6.47,6.30e2,2.53,3.79e3,8.99e3,6.75e2},{0.3750,1.74,3.15,2.78,6.37,7.06e2,2.53,3.67e3,1.02e4,7.35e2},{0.3618,1.78,3.14,2.78,6.33,7.23e2,2.53,3.37e3,1.06e4,7.44e2},{0.3524,1.80,3.13,2.79,6.29,7.57e2,2.53,3.28e3,1.11e4,7.69e2},{0.3453,1.82,3.12,2.79,6.26,7.94e2,2.53,3.25e3,1.17e4,7.93e2},{0.3399,1.84,3.11,2.79,6.23,8.17e2,2.53,3.20e3,1.2e4,8.25e2}};
const double betavalues[8][11] = {{0.7029,1.05e3,6.33,3.71e3,4.11e3,5.15e2,2.54,1.5,2.87,3.27,7.09},{0.9054,1.38e3,6.33,3.8e3,7.05e3,6.42e3,4.40,1.5,2.43,5.18,9.34},{1.018,1.55e3,6.33,3.8e3,8.75e3,6.9e2,3.77,1.5,1.46,4.34,8.65},{1.092,1.64e3,6.33,3.74e3,9.73e3,7.07e2,3.43,1.5,1.06,3.92,8.27},{1.146,1.71e3,6.33,3.72e3,1.07e4,7.31e2,3.2,1.5,0.848,3.66,8.02},{1.186,1.74e3,6.33,3.64e3,1.11e4,7.35e2,3.05,1.5,0.718,3.38,7.83},{1.218,1.73e3,6.33,3.53e3,1.13e4,7.29e2,2.92,1.5,0.629,3.33,7.68},{1.244,1.79e3,6.33,3.57e3,1.17e4,7.34e2,2.82,1.5, 0.565,3.21,7.55}};
const double kappavalues[8][11] = {{3.203,6.18,4.66,1.93,2.31,5.35,4.01,2.5,0.661,0.931,2.5},{4.931,9.3,3.96,1.89,3.78,7.78,2.46,2.5,0.156,0.398,1.71},{6.115,11,3.72,1.66,4.76,8.88,1.13,2.5,0.0442,0.175,1.05},{6.995,9.14,3.6,1.31,4.63,8.8,0.628,2.5,0.018,0.101,0.775},{7.68,8.6,3.53,1.12,4.62,8.8,0.418,2.5,0.00963,0.0702,0.646},{8.231,8.57,3.49,1.04,4.83,8.96,0.319,2.5,0.00625,0.0551,0.578},{8.685,8.84,3.49,1.02,5.19,9.24,0.268,2.5,0.00461,0.0465,0.539},{9.067,7.93,3.43,0.875,4.74,8.84,0.238,2.5,0.00371,0.041,0.515}};


//The first transport coefficients are the RESISTIVITY parallel to the B field, perpendicular to the plane of the
//field and transport 'force' and the wedge is the resistivity perpendicular to the previous directions

double alpha_para(double chi, int Z){

    double alpha0 = alphavalues[Z-1][0];

    return alpha0/brag;
}

double alpha_perp(double chi, int Z){


    double alpha0p = alphavalues[Z-1][1];
    double alpha1p = alphavalues[Z-1][2];
    double a0p = alphavalues[Z-1][3];
    double a1p = alphavalues[Z-1][4];


    return (1-(alpha1p*chi+alpha0p)/(chi*chi+a1p*chi+a0p))/brag;

}

double alpha_wedg(double chi, int Z){

    double alpha0pp = alphavalues[Z-1][5];
    double alpha1pp = alphavalues[Z-1][6];
    double a0pp = alphavalues[Z-1][7];
    double a1pp = alphavalues[Z-1][8];
    double a2pp = alphavalues[Z-1][9];

    return -chi*(alpha1pp*chi+alpha0pp)/pow((chi*chi*chi+a2pp*chi*chi+a1pp*chi+a0pp),double(8.0/9.0))/brag;
}

//The following are the THERMOELECTRIC components

double beta_para(double chi, int Z){

    double beta0 = betavalues[Z-1][0];

    return beta0;
}

double beta_perp(double chi, int Z){

    double beta0p = betavalues[Z-1][1];
    double beta1p = betavalues[Z-1][2];
    double b0p = betavalues[Z-1][3];
    double b1p = betavalues[Z-1][4];
    double b2p = betavalues[Z-1][5];

    return (beta1p*chi+beta0p)/pow((chi*chi*chi+b2p*chi*chi+b1p*chi+b0p),double(8.0/9.0));
}

double beta_wedg(double chi, int Z){

    double beta0pp = betavalues[Z-1][6];
    double beta1pp = betavalues[Z-1][7];
    double b0pp = betavalues[Z-1][8];
    double b1pp = betavalues[Z-1][9];
    double b2pp = betavalues[Z-1][10];

    return chi*(beta1pp*chi+beta0pp)/(chi*chi*chi+b2pp*chi*chi+b1pp*chi+b0pp);
}

//The following are the THERMAL CONDUCTIVITY components

double kappa_para(double chi, int Z){

    double gamma0 = kappavalues[Z-1][0];

    return brag*gamma0;
}


double kappa_perp(double chi, int Z){

    double gamma0p = kappavalues[Z-1][1];
    double gamma1p = kappavalues[Z-1][2];
    double c0p = kappavalues[Z-1][3];
    double c1p = kappavalues[Z-1][4];
    double c2p = kappavalues[Z-1][5];

    return brag*(gamma1p*chi+gamma0p)/(chi*chi*chi+c2p*chi*chi+c1p*chi+c0p);
}


double kappa_wedg(double chi, int Z){

    double gamma0pp = kappavalues[Z-1][6];
    double gamma1pp = kappavalues[Z-1][7];
    double c0pp = kappavalues[Z-1][8];
    double c1pp = kappavalues[Z-1][9];
    double c2pp = kappavalues[Z-1][10];

    return brag*chi*(gamma1pp*chi+gamma0pp)/(chi*chi*chi+c2pp*chi*chi+c1pp*chi+c0pp);
}

//#####################################################
//
//     HALL PARAMETER: Function for evaluating the hall
//     parameter at mesh position i,j,k, this will be used to
//     calculate the transport coefficients
//
//#####################################################

double hall(int Z, double *vn, int position){

    return brag*sqrt(vn[position+8]*vn[position+8]+vn[position+9]*vn[position+9]+vn[position+10]*vn[position+10])*(pow(vn[position+4],1.5)/vn[position]);

}

//###################################################
//
//    MAGNETIC FIELD DIRECTION: calculates the ratio of
//    a compnent to the norm of the Bfield
//
//###################################################

void setbdir(int position, double *vn, double *b){

  double Bmag=sqrt(vn[position+8]*vn[position+8]+vn[position+9]*vn[position+9]+vn[position+10]*vn[position+10]);
  if (Bmag<0.001) {
    b[0]=0;
    b[1]=0;
    b[2]=0;     //This is a cutoff so to avoid causingproblems by doing 0/0
  }
  else{
    b[0]=vn[position+8]/Bmag;
    b[1]=vn[position+9]/Bmag;
    b[2]=vn[position+10]/Bmag;
  }

}

//#########################################################################################
//
//                   CARTESIAN COMPONENTS OF TRANSPORT COEFFICIENTS
//    The following are functions to set the 9 cartesian compnents of the transport
//    coefficients. It is easier in 3D to code up using cartesian components xx, xy, xz etc
//    which is what these functions compute using the magnetic
//    field and the functions above.
//
//#########################################################################################

//coefficient tensor
void setcoeff(int Z, double chi, double *out,double *b,int q){

  double a1,a2,a3;
  if (q==1) {
    a1=alpha_para(chi,Z);
    a2=alpha_perp(chi, Z);
    a3=alpha_wedg(chi,Z);
  }
  else if (q==2) {
    a1=beta_para(chi,Z);
    a2=beta_perp(chi, Z);
    a3=beta_wedg(chi,Z);
  }
  else{
    a1=kappa_para(chi,Z);
    a2=kappa_perp(chi, Z);
    a3=kappa_wedg(chi,Z);
  }

  out[0]=a2+(a1-a2)*b[0]*b[0];
  out[1]=-a3*b[2]+(a1-a2)*b[0]*b[1];
  out[2]=a3*b[1]+(a1-a2)*b[0]*b[2];
  out[3]=a3*b[2]+(a1-a2)*b[0]*b[1];
  out[4]=a2+(a1-a2)*b[1]*b[1];
  out[5]=-a3*b[0]+(a1-a2)*b[1]*b[2];
  out[6]=-a3*b[1]+(a1-a2)*b[0]*b[2];
  out[7]=a3*b[0]+(a1-a2)*b[1]*b[2];
  out[8]=a2+(a1-a2)*b[2]*b[2];

}
