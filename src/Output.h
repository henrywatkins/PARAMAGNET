//####################################################
//
//                  OUTPUT HEADER
//
//     This is just the header that is paired with
//     the Output.cpp data file which contains the
//              functions for output
//
//####################################################


#ifndef Output_h
#define Output_h

void ASCIIOutput(std::ofstream& myfile, int var, double *y );
void ASCIIOutputvec(std::ofstream& myfile, int var, double *y );
void ASCIILaser(std::ofstream& myfile, std::complex<double> *z);
void ASCIILasercorr(std::ofstream& myfile, std::complex<double> *z);
void VTKOutput(std::ofstream& parafile, double *y,std::complex<double> *z,std::complex<double> *z2);

#endif /* Output_hpp */
