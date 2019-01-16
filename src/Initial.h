//####################################################
//
//                   INITIAL HEADER
//
//     This is just the header that is paired with
//     the Initial.cpp file which contains the initial
//     condition setter and laser profile setter
//
//####################################################


#ifndef LaserProfile_h
#define LaserProfile_h

void InitialCondition(double *initvalue);
void InitialLaser(std::complex<double> *initvalue,std::complex<double> *initvalue2,int t);

#endif
