//###################################################################
//
//                              OUTPUT
//
//       This file contains functions for producing the output, it
//       produces simple ASCII text file type and legacy VTK formats
//
//
//###################################################################
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
#include "Output.h"
#include "Global.h"
#include <cmath>

using namespace std;

//ASCIIOutput

void ASCIIOutput(std::ofstream& myfile, int var, double *y){

    for (int index=0; index<N*vars; index=index+vars) {
        myfile << y[index+var-1] <<" ";
    }
    myfile<<endl;
}

void ASCIILaser(std::ofstream& myfile, std::complex<double> *z){

  for (int index=0; index<N; index++) {
      myfile << norm(z[index])<<" ";
  }
  myfile<<endl;
}

void ASCIILasercorr(std::ofstream& myfile, std::complex<double> *z){

  for (int index=0; index<N; index++) {
      myfile << z[index]<<" ";
  }
  myfile<<endl;
}

void ASCIIOutputvec(std::ofstream& myfile, int var, double *y){

    for (int index=0; index<N*vars; index=index+vars) {
        myfile << y[index+var-1] <<" "<<endl;
        myfile << y[index+var] <<" "<<endl;
        myfile << y[index+var+1] <<" "<<endl;
    }
    myfile<<endl;
}

//VTK output
void VTKOutput(std::ofstream& parafile, double *y,std::complex<double> *z,std::complex<double> *z2){

    int p;

    parafile <<"# vtk DataFile Version 3.0"<<endl;
    parafile <<"3D MAGNET output"<<endl;
    parafile <<"ASCII"<<endl;
    parafile <<"DATASET STRUCTURED_GRID"<<endl;
    parafile <<"DIMENSIONS "<<nx-2<<" "<<ny-2<<" "<<nz-2<<endl;
    parafile <<"POINTS "<< (nx-2)*(ny-2)*(nz-2)<<" double"<<endl;
    for (int k=0; k<nz-2; k++) {
      for (int j=0; j<ny-2; j++){
          for (int i=0; i<nx-2; i++) {
              parafile <<i*dx<<" "<<j*dy<<" "<<k*dz<<endl;
          }
      }
    }
    parafile <<"POINT_DATA "<<(nx-2)*(ny-2)*(nz-2)<<endl;

    //Density
    if (nout==1) {
        parafile <<"SCALARS Density double "<<endl;
        parafile <<"LOOKUP_TABLE default"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p]<<endl;
                }
            }
        }

    }
    //Velocity
    if (Vout==1) {
        parafile <<"VECTORS Velocity double"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p+1]<<" "<<y[p+2]<<" "<<y[p+3]<<endl;
                }
            }
        }
    }

    //Temperature
    if (Tout==1) {
        parafile <<"SCALARS Temperature double "<<endl;
        parafile <<"LOOKUP_TABLE default"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p+4]<<endl;
                }
            }
        }
    }

    //Electric Field
    if (Eout==1) {
        parafile <<"VECTORS ElectricField double"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p+5]<<" "<<y[p+6]<<" "<<y[p+7]<<endl;
                }
            }
        }
    }
    //Magnetic Field
    if (Bout==1) {
        parafile <<"VECTORS MagneticField double"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p+8]<<" "<<y[p+9]<<" "<<y[p+10]<<endl;
                }
            }
        }
    }
    //Current
    if (jout==1) {
        parafile <<"VECTORS Current double"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p+11]<<" "<<y[p+12]<<" "<<y[p+13]<<endl;
                }
            }
        }
    }
    //Heat Flow
    if (qout==1) {
        parafile <<"VECTORS HeatFlow double"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=vars*(i+nx*j+nx*ny*k);
                  parafile <<y[p+14]<<" "<<y[p+15]<<" "<<y[p+16]<<endl;
                }
            }
        }
    }

    if (lout==1) {
        parafile <<"SCALARS LaserField double"<<endl;
        parafile <<"LOOKUP_TABLE default"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=(i+nx*j+nx*ny*k);
                  parafile <<norm(z[p])<<endl;
                }
            }
        }
    }
    if (l2out==1) {
        parafile <<"SCALARS LaserField2 double"<<endl;
        parafile <<"LOOKUP_TABLE default"<<endl;
        for (int k=1; k<nz-1; k++) {
            for (int j=1; j<ny-1; j++){
                for (int i=1; i<nx-1; i++) {
                  p=(i+nx*j+nx*ny*k);
                  parafile <<norm(z2[p])<<endl;
                }
            }
        }
    }

    parafile.close();

}
