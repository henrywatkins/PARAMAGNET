//###################################################################
//
//                            INITIAL DATA
//
//        This file contains the functions for setting the initial
//        condition and the laser beams incident in the plasma.
//        It also contains the data associated with the laser heating
//        such as wavelength, intensity and beam profile.
//
//###################################################################

#include <complex>
#include <math.h>
#include "Global.h"
#include "Initial.h"
#include <cstdlib>
#include <fstream>
#include <sstream>           //String functions
#include <string>

using namespace std;

//if the inital value of a variable is not defined it will automattically be set to zero
void InitialCondition(double *initvalue){

   int i,j,k;
   int l;                             //row indexing index
   int cl;                            //C index of the vector, since C goes from 0 to N-1 rather than 1 to N
   int nin,Vxin,Vyin,Vzin,Tin,Exin,Eyin,Ezin,Bxin,Byin,Bzin,jxin,jyin,jzin,qxin,qyin,qzin, liin,lrin;


   stringstream ss;
   double array[nx][ny][nz];
   string filename="./profile.txt";
   ifstream file(filename.c_str());
   if (file.is_open()) {
     for (size_t k = 0; k < nz; k++) {
       for (size_t j = 0; j < ny; j++) {
         for (size_t i = 0; i < nx; i++) {
           file >> array[i][j][k];
         }
       }
     }
   }
   file.close();

   for (k=0; k<nz; k++) {
       for (j=0; j<ny; j++) {
           for (i=0; i<nx; i++) {

               //Set index value for row indexed scheme
               l=vars*(i+nx*j+nx*ny*k);

               //set index value for each variable
               nin=l;
               Vxin=l+1;
               Vyin=l+2;
               Vzin=l+3;
               Tin=l+4;
               Exin=l+5;
               Eyin=l+6;
               Ezin=l+7;
               Bxin=l+8;
               Byin=l+9;
               Bzin=l+10;
               jxin=l+11;
               jyin=l+12;
               jzin=l+13;
               qxin=l+14;
               qyin=l+15;
               qzin=l+16;

               //#########################################
               //
               //     INITIAL CONDITION SETTING STAGE
               //  To set an initial condition this section
               //  and only this section should be modified,
               //  each variable can be set independently,
               //
               //#########################################

               initvalue[nin]=1.0;//+(i*dx*0.93/Lx+k*dz*0.36/Lz);
               initvalue[Tin]=1.0+0.01*array[i][j][k];
               //initvalue[Vxin]=-0.1*sin(2*3.14*xpara)*sn;
               //initvalue[Vyin]=0.1*sin(2*3.14*xpara)*cs;
               //initvalue[Vzin]=0.1*cos(2*3.14*xpara);
               //initvalue[Bxin]=7.5*tei0*(1.75e11)*0.36;
               //initvalue[Byin]=5*tei0*(1.75e11)*(1+k*dz/Lz);///sqrt(2);
               initvalue[Bzin]=5*tei0*(1.75e11);
           }
       }
   }

}

void InitialLaser(std::complex<double> *initvalue, std::complex<double> *initvalue2, int t){

 //Gaussian Beam Variables
 double Rlength, phi, waist,waist0,Rz, sep;

 //waist0=77e-6;
 //Rlength=3.14*waist0*waist0/wl;
 //Rz=-0.5*zDom*(1+pow(Rlength/0.5/zDom,2.0));
 //phi=atan(-0.5*zDom/Rlength);
 //waist=500e-6/l0;//waist0*sqrt((1+pow(0.5*zDom/Rlength,2.0)));
 waist=400e-6/2.35/l0;

 int l;
 double re,im;
 for (int j = 1; j < ny-1; j++) {
   for (int i = 1; i < nx-1; i++) {

     l=i+nx*j+nx*ny;
     initvalue[l]=complex<double>(exp(-0.5*pow((i*dx-nx*dx/2)/waist,2.0)-0.5*pow((j*dx-ny*dy/2)/waist,2.0)),0.0);
     initvalue2[l]=complex<double>(0.0,0.0);

                             //3D//waist0*exp(-pow((i-nx/2)*dx*l0/waist,2.0)-pow((j-ny/2)*dy*l0/waist,2.0))*cos(0.5*wn*(pow((i-nx/2)*dx,2.0)+pow((j-ny/2)*dy,2.0))*l0/Rz-phi)/waist,waist0*exp(-pow((i-nx/2)*dx*l0/waist,2.0)-pow((j-ny/2)*dy*l0/waist,2.0))*sin(0.5*wn*(pow((i-nx/2)*dx,2.0)+pow((j-ny/2)*dy,2.0))*l0/Rz-phi)/waist
                             //2D//waist0*exp(-pow((i-nx/2)*dx*l0/waist,2.0))*cos(0.5*wn*(pow((i-nx/2)*dx,2.0))*l0/Rz-phi)/waist,waist0*exp(-pow((i-nx/2)*dx*l0/waist,2.0))*sin(0.5*wn*(pow((i-nx/2)*dx,2.0))*l0/Rz-phi)/waist
                             //array[i][j],0.0
   }                         //(1+0.01*im)*exp(-pow((i-nx/2)*dx/waist,4.0)),0.01*re*exp(-pow((i-nx/2)*dx/waist,4.0))
 }
}
