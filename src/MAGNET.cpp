//#########################################################################################################################
//
//                 ██████╗  █████╗ ██████╗  █████╗ ███╗   ███╗ █████╗  ██████╗ ███╗   ██╗███████╗████████╗
//                 ██╔══██╗██╔══██╗██╔══██╗██╔══██╗████╗ ████║██╔══██╗██╔════╝ ████╗  ██║██╔════╝╚══██╔══╝
//                 ██████╔╝███████║██████╔╝███████║██╔████╔██║███████║██║  ███╗██╔██╗ ██║█████╗     ██║
//                 ██╔═══╝ ██╔══██║██╔══██╗██╔══██║██║╚██╔╝██║██╔══██║██║   ██║██║╚██╗██║██╔══╝     ██║
//                 ██║     ██║  ██║██║  ██║██║  ██║██║ ╚═╝ ██║██║  ██║╚██████╔╝██║ ╚████║███████╗   ██║
//                 ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝   ╚═╝
//
//
//                                (PARAllel MAGnetized Newton-method code for Electron Transport)
//
//                                       MAGNETIZED COLLISIONAL TRANSPORT CODE
//
//                         Author: Henry Watkins, The Blackett Laboratory, Imperial College London
//                                                   Version: 6.0
//
//           Description:
//           The code is a 3D MHD plasma code with full electron collisional transport. It uses an  implicit
//           finite difference scheme for the plasma transport equations and solves using finite differences with Newton's method.
//
//           Compiling and running:
//           This code relies quite heavily on the PETSc libraries and does not compile without the particular
//           libraries of PETSc version 3.6.3 or above. To compile first see the user guide to ensure the required
//           libraries and files are on the computer then use the command 'make MAGNET' in the terminal (working within
//           the PARAMAGNET directory). To run in serial simply use './MAGNET' however if the parallelized version is being
//           used run with 'mpiexec -n (number of cores) ./MAGNET'
//
//
//           Output:
//           This code produces two kinds of output; it produces a .vtk database that can be read
//           into Paraview or VisIt (with some specification of the filters). The second type are ASCII text files
//           for each variable that can be read into matlab or matplotlib. The output format can be read using
//           several matlab functions inluded in the source code directory.
//
//           Method:
//           This code solves the finite differenced set of PDEs describing collisional MHD with full electron transport.
//           In this version all variables are treated implicitly using a backward euler time differencing and the space
//           derivatives are all central space differences. The result is the matrix equation:
//
//                                                         A(x)x=B(x)b,
//
//           where x is the sytem at timestep n+1 and b is the system at timestep n, the matricies are both functions of x.
//           This nonlinear matrix equation is solved using the multivariable Newton/Raphson method within which a Krylov
//           subspace method is used to solve the linear problem.
//
//
//
//##########################################################################################################################


//################################################################################################
//
//      SECTION 1: Headers and function delarations
//
//################################################################################################

#include <complex>
#include <petscsnes.h>       //inludes all the petsc solvers and all other petsc utilities
#include <iostream>          //Console I/O
#include <fstream>           //File I/O
#include <math.h>            //Math functions
#include <petsctime.h>       //inludes a timer for the program
#include <sstream>           //String functions
#include <string>            //More sting functions important for the VTK output file routines
#include "Global.h"          //inludes all the parameters specific to the problem set by the user in the Global.cpp file
#include "TransportCoeff.h"  //inludes the functions for evaluation the transport coefficients
#include "Output.h"          //Inludes functions for the output
#include "Initial.h"         //Initial condition setter
#include "omp.h"


using namespace std;

struct datastruct{
    Vec        b;             //The RHS vector b as well as the laser vector
    complex<double> *laser;
    complex<double> *laser2;
};

extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);     //delaration of the function evaluation function of section 4
extern PetscErrorCode FormInitialGuess(SNES, Vec, void*);   //delaration of the inital guess of the vector x in section 6
PetscErrorCode PARAXIAL(datastruct &ctx,complex<double>* );
void zTriSolver(complex<double> , complex<double>*, complex<double> , complex<double> *,int );
void SM(complex<double> , complex<double>*, complex<double> , complex<double> *,int );

//################################################################################################
//
//      SECTION 2: Begin main program and set the solution type
//
//################################################################################################

int main(int argc,char **argv)
{

    //###################################################
    //
    //      SECTION 2.1: Data/Structure initialisation
    //
    //###################################################


    PetscLogDouble     t1,t2,t3,t4,t5,elapsed_time;   //Variables for timing the program
    SNES               snes;                 //nonlinear solver context
    KSP                ksp;                  //KSP context-allows you to hard code which KSP method
    PC                 pc;                   //PC context-allows you to hard code which PC method
    Vec                x,f,b;                //solution, function vector , Ax=Bb
    Mat                J;                    //Jacobian matrix
    PetscInt           its,its2;                  //Variable to store the number of iterations before the solutions converges within the tolerance
    datastruct         user;                 //The parallel data held in the structure above

    user.laser=new complex<double>[N];
    user.laser2=new complex<double>[N];

    PetscInitialize(&argc,&argv,(char*)0,PETSC_NULL);   //Initialize the program, required for all petsc programs

    PetscTime(&t1);                      //First time point for timing the program

    cout<<"alpha= "<<alpha<<" beta= "<<beta<<" gamma= "<<gammad<<" IBC="<<IBc<<" Pc="<<Pc<<endl;

    SNESCreate(PETSC_COMM_WORLD,&snes);  //Create nonlinear solve

    //Create vectors, this must be done for any structure used in the code

    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,N*vars);  //the size of the solution vector will be nx*ny*nz*vars
    VecSetFromOptions(x);                //Required, also allows certain options to be called from the command line
    VecDuplicate(x,&b);
    user.b = b;               //likewise for the previous timestep vector
    VecDuplicate(x,&f);
    VecSet(f,0);


    //###################################################
    //
    //      SECTION 2.2: Initial Conditions, set b, set laser
    //
    //###################################################


    //Makes a vector with length N*vars where the elements are the indicies, requied for pulling values out of the solution vectors
    PetscInt *ix=new PetscInt[N*vars];
    for (int index=0; index<N*vars; index++) {
        ix[index]=index;
    }

    double *initvalue=new double[N*vars];

    //PetscScalar initvalue[N*vars]; //initialise vector to hold values

    InitialCondition(initvalue);  //Calls the function from the Initial.cpp files that fills the initial state

    VecSetValues(b,N*vars,ix,initvalue,INSERT_VALUES);  //Sets the values of the Petsc vector b to the initial state

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    PetscPrintf(PETSC_COMM_WORLD,"Initial Conditions Set \n"); //use PetscPrintf instead of cout because in parallel cout will print as many times as there are processes

    //InitialLaser(user.laser);

    delete[] initvalue;
    delete[] ix;


    //###################################################
    //
    //      SECTION 2.3: Solution method of the nonlinear
    //    matrix. This allows you to choose the SNES
    //    nonlinear method and Krylov subspace and
    //    preconditioner for solving the equation. The
    //    SNES method list is on p106 of Petsc manual,
    //    the KSP on p77 and preconditoners on p80.
    //    The SNES method is set to the newton method
    //    and should not be changed.
    //
    //###################################################


    //The current options set for the PC, KSP methods have been found to give the best performance on the tests I have run. They are however optional.

    SNESGetKSP(snes, &ksp);
    SNESSetTolerances(snes,1e-9,1e-9,PETSC_DEFAULT,PETSC_DEFAULT,500);
    KSPSetType(ksp, KSPGMRES);              //Choose the KSP linear method being used
    KSPSetTolerances(ksp,1e-10,1e-10,PETSC_DEFAULT,2000);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCNONE);
    SNESSetType(snes,SNESNEWTONLS);         //Set the nonlinear method being ued to solve the problem
    SNESSetFromOptions(snes);

    //################################################################################################
    //
    //      SECTION 3: Solution
    //
    //################################################################################################


    //###################################################
    //
    //      SECTION 3.1: Set up file output
    //
    //###################################################

    //initialises a vector of length N*vars that contains the data of the solution once it has been pulled out of the solution vector

    PetscScalar *y;

    //pulls values out of b and puts them into y so that it can be outputted to a file, same with the laser vector z
    VecGetArray(user.b,&y);

    //initialise time point, this is the last timestep where data was outputted
    int tlast=0;


    //###################################################
    //
    //      SECTION 3.1.1: Set up ASCII format output
    //      using function from Output.cpp file
    //
    //###################################################

    //The if statements come from the variable output choices made in the Global.cpp file

    //Density output
    ofstream myfile1;
    myfile1.open("./ASCIIOutput/noutput.txt");

    //Velocity output
    ofstream myfile2;
    myfile2.open("./ASCIIOutput/Voutput.txt");

    //Temperature output
    ofstream myfile3;
    myfile3.open("./ASCIIOutput/Toutput.txt");

    //Electric field output
    ofstream myfile4;
    myfile4.open("./ASCIIOutput/Eoutput.txt");

    //Magnetic field output
    ofstream myfile5;
    myfile5.open("./ASCIIOutput/Boutput.txt");

    //Current output
    ofstream myfile6;
    myfile6.open("./ASCIIOutput/joutput.txt");

    //Heat flow output
    ofstream myfile7;
    myfile7.open("./ASCIIOutput/qoutput.txt");

    //laser output
    ofstream myfile8;
    myfile8.open("./ASCIIOutput/laserintensity.txt");

    ofstream myfile9;
    myfile9.open("./ASCIIOutput/laserintensity2.txt");



    //###################################################
    //
    //      SECTION 3.1.2: Set up VTK format file output
    //      using output function from Ouput.cpp, then
    //      lose the pull out with VecRestoreArray
    //
    //###################################################

    ofstream parafile;
    parafile.open("./VTKOutput/VTKoutput.0.vtk");

    VTKOutput(parafile,y,user.laser,user.laser2);

    VecRestoreArray(user.b,&y);


    //###################################################
    //
    //      SECTION 3.2: Begin the time integration and solve
    //
    //###################################################

    for (int t=0; t<m; t++) {


      //##################################################
      //
      //    SECTION 3.4: Laser solver
      //
      //##################################################
      InitialLaser(user.laser,user.laser2,t);
      PARAXIAL(user,user.laser);
      //PARAXIAL(user,user.laser2);


      //###################################################
      //
      //    SECTION 3.5: Form the residual vector f and
      //    solve using the matrix free method from the
      //    PETSc SNES library
      //
      //###################################################

      //Set function evaluation routine and vector. Has to be done for each timestep

      SNESSetFunction(snes,f,FormFunction,&user);

      //Set Jacobian matrix data structure and Jacobian evaluation routine, done for each timestep
      MatCreateSNESMF(snes,&J);

      SNESSetJacobian(snes,J,J,MatMFFDComputeJacobian,&user);

      //Set the inital guess of x, this should be done before calling SNESSolve

      SNESSetComputeInitialGuess(snes,FormInitialGuess,&user);

      //Solve the system for x
      SNESSolve(snes,NULL,x);

      PetscPrintf(PETSC_COMM_WORLD,"Timestep Solved\n");

      //Find the number of iterations required to reach convergence

      SNESGetIterationNumber(snes,&its);
      PetscPrintf(PETSC_COMM_WORLD,"Number of nonlinear iterations for this timestep = %D\n",its);

     //###################################################
     //
     //    SECTION 3.6: Set the solution x to be the b for
     //    the next timestep and check convergence, this
     //    is a small but very important step since this timesteps the problem
     //
     //###################################################

      VecAYPX(b,0,x);                 //sets old vector to new values

      PetscPrintf(PETSC_COMM_WORLD,"Completed Timestep %d/%d\n",t+1,m);

      MatDestroy(&J);

      //###################################################
      //
      //      SECTION 3.3: Output data for this timestep
      //
      //###################################################


      //choose only a sample of the timesteps to output, otherwise data requirement can get very large
      if (t==tlast) {

          //Pull the data out of the vector x and stick it in y
          VecGetArray(user.b,&y);

          PetscPrintf(PETSC_COMM_WORLD,"Outputting Data\n");

          //###################################################
          //
          //      SECTION 3.3.1: VTK Output data for this timestep
          //
          //###################################################
          //set filename

          ostringstream filename;
          filename<<"./VTKOutput/VTKoutput."<<t<<".vtk";

          //open file and output

          ofstream parafile;
          parafile.open(filename.str().c_str());


          VTKOutput(parafile, y,user.laser,user.laser2);

          //###################################################
          //
          //      SECTION 3.3.2: ASCII output data for this timestep
          //
          //###################################################

          //if statements choose whether to output that particular variable, defined in the output section of Global.cpp

          //Density output
          if(nout==1 ){ASCIIOutput(myfile1,1,y);}

          //Velocity in x direction output
          if(Vout==1 ){ASCIIOutputvec(myfile2,3,y);}

          //Temperature output
          if(Tout==1 ){ASCIIOutput(myfile3,5,y);}

          //Electric field in x direction output
          if(Eout==1 ){ASCIIOutputvec(myfile4,7,y);}

          //Magnetic field in x direction output
          if(Bout==1 ){ASCIIOutputvec(myfile5,10,y);}

          //Current in x direction output
          if(jout==1 ){ASCIIOutputvec(myfile6,13,y);}

          //Heat flow output
          if(qout==1 ){ASCIIOutputvec(myfile7,16,y);}

          //laser output
          if(lout==1 ){ASCIILaser(myfile8,user.laser);}

          if(l2out==1 ){ASCIILaser(myfile9,user.laser2);}

          VecRestoreArray(x,&y);

          //update the last timepoint where data was output
          tlast+=tfrac;

      }
      SNESGetIterationNumber(snes,&its);
      if (its==0) {

          PetscPrintf(PETSC_COMM_WORLD, "The system appears to have stopped changing, this could be a physical steady state or more likely numerical collapse. Evaluation Terminating\n");

          return 0;
      }


}

    //###################################################
    //
    //    SECTION 3.7: Finish up the program. lose the
    //    ASCII output files, finish timer and free work space.  All
    //    PETSc objects should be destroyed when they
    //    are no longer needed.
    //
    //###################################################

    //lose output files

    myfile1.close();
    myfile2.close();
    myfile3.close();
    myfile4.close();
    myfile5.close();
    myfile6.close();
    myfile7.close();
    myfile8.close();
    myfile9.close();
    //Time measurement section for the program duration

    PetscTime(&t2); //Second time point
    elapsed_time = t2-t1;  //elapsed time
    PetscPrintf(PETSC_COMM_WORLD,"Program Duration: %f seconds\n",elapsed_time); //outputs time duration of program


    //destory petsc objects
    delete[] user.laser;
    delete[] user.laser2;
    VecDestroy(&x);
    VecDestroy(&f);
    VecDestroy(&b);
    SNESDestroy(&snes);
    PetscFinalize();
    return 0;
}


//################################################################################################
//
//      SECTION 4: Function evaluation
//
//      FormFunction - Evaluates nonlinear function, f(x)= A(x)x-B(x)b
//
//      Input Parameters:
//      .  snes - the SNES context
//      .  x    - solution vector
//      .  ctx - allows the RHS b to be an input
//
//      Output Parameter:
//      .  f - function vector
//
//################################################################################################

PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void *ctx){

    //###################################################
    //
    //      SECTION 4.1: Data initialisation
    //
    //###################################################
    Vec ref;
    //const PetscScalar *x2;             //solution vector x
    PetscScalar *xx;
    PetscScalar *bb;             //previous timestep solution vector b
    PetscScalar *ff;                   //nonlinear function vector f
                               //row indexing index
    int d1x,d2x,d1y,d2y,d1z,d2z;       //the difference values that allow easy application of boundary conditions

    datastruct       *user = (datastruct*) ctx;
    Vec              b = user->b;
    complex<double>  *laser=user->laser;
    complex<double>  *laser2=user->laser2;

    //Get pointers to vector data.
    VecDuplicate(x, &ref);
    VecCopy(x,ref);
    VecGetArray(ref,&xx);
    VecGetArray(b,&bb);
    VecGetArray(f,&ff);

    d1x=vars;
    d2x=-vars;
    d1y=vars*nx;
    d2y=-vars*nx;
    d1z=vars*nx*ny;
    d2z=-vars*nx*ny;

    //###################################
    //
    //    initialise variables so they
    //           are thread safe
    //
    //###################################

    #pragma omp parallel
    {
      int l;
      double bdir[3];
      double resistivity[9];
      double thermoelectric[9];
      double conductivity[9];
      int nin,Vxin,Vyin,Vzin,Tin,Exin,Eyin,Ezin,Bxin,Byin,Bzin,jxin,jyin,jzin,qxin,qyin,qzin;      //total index, the result of both the row indexing index and the variable index combined
      double chi,bierx,biery,bierz, qlimx,qlimy,qlimz,qmag,qfs, lasnorm,lasnormx1,lasnormx2,lasnormy1,lasnormy2,lasnormz1,lasnormz2;

    //##################################################
    //
    //      SECTION 4.2: Boundary Conditions
    //
    //##################################################
    #pragma omp for
    for (int k = 1; k < nz-1; k++) {
      for (int j = 1; j < ny-1; j++) {
        for (int a = 0; a < vars; a++) {
          if (xbcon==1) {
            //Relfective
            xx[a+vars*(nx*j+nx*ny*k)]=xx[a+vars*(1+nx*j+nx*ny*k)];
            xx[a+vars*(nx-1+nx*j+nx*ny*k)]=xx[a+vars*(nx-2+nx*j+nx*ny*k)];
            bb[a+vars*(nx*j+nx*ny*k)]=bb[a+vars*(1+nx*j+nx*ny*k)];
            bb[a+vars*(nx-1+nx*j+nx*ny*k)]=bb[a+vars*(nx-2+nx*j+nx*ny*k)];
          }
          else if (xbcon==3) {
            //outflow
            xx[a+vars*(nx*j+nx*ny*k)]=2*xx[a+vars*(1+nx*j+nx*ny*k)]-xx[a+vars*(2+nx*j+nx*ny*k)];
            xx[a+vars*(nx-1+nx*j+nx*ny*k)]=2*xx[a+vars*(nx-2+nx*j+nx*ny*k)]-xx[a+vars*(nx-3+nx*j+nx*ny*k)];
            bb[a+vars*(nx*j+nx*ny*k)]=2*bb[a+vars*(1+nx*j+nx*ny*k)]-bb[a+vars*(2+nx*j+nx*ny*k)];
            bb[a+vars*(nx-1+nx*j+nx*ny*k)]=2*bb[a+vars*(nx-2+nx*j+nx*ny*k)]-bb[a+vars*(nx-3+nx*j+nx*ny*k)];
          }
          else{
            //Periodic
            xx[a+vars*(nx*j+nx*ny*k)]=xx[a+vars*(nx-2+nx*j+nx*ny*k)];
            xx[a+vars*(nx-1+nx*j+nx*ny*k)]=xx[a+vars*(1+nx*j+nx*ny*k)];
            bb[a+vars*(nx*j+nx*ny*k)]=bb[a+vars*(nx-2+nx*j+nx*ny*k)];
            bb[a+vars*(nx-1+nx*j+nx*ny*k)]=bb[a+vars*(1+nx*j+nx*ny*k)];
          }
        }
      }
    }


    //yBoundaries
    #pragma omp for
    for (int k = 1; k < nz-1; k++) {
      for (int i = 1; i < nx-1; i++) {
        for (int a = 0; a < vars; a++) {
          if (ybcon==1) {
            //Relfective
            xx[a+vars*(i+nx*ny*k)]=xx[a+vars*(i+nx+nx*ny*k)];
            xx[a+vars*(i+nx*(ny-1)+nx*ny*k)]=xx[a+vars*(i+nx*(ny-2)+nx*ny*k)];
            bb[a+vars*(i+nx*ny*k)]=bb[a+vars*(i+nx+nx*ny*k)];
            bb[a+vars*(i+nx*(ny-1)+nx*ny*k)]=bb[a+vars*(i+nx*(ny-2)+nx*ny*k)];
          }
          else if (ybcon==3) {
            //outflow
            xx[a+vars*(i+nx*ny*k)]=2*xx[a+vars*(i+nx+nx*ny*k)]-xx[a+vars*(i+nx*2+nx*ny*k)];
            xx[a+vars*(i+nx*(ny-1)+nx*ny*k)]=2*xx[a+vars*(i+nx*(ny-2)+nx*ny*k)]-xx[a+vars*(i+nx*(ny-3)+nx*ny*k)];
            bb[a+vars*(i+nx*ny*k)]=2*bb[a+vars*(i+nx+nx*ny*k)]-bb[a+vars*(i+nx*2+nx*ny*k)];
            bb[a+vars*(i+nx*(ny-1)+nx*ny*k)]=2*bb[a+vars*(i+nx*(ny-2)+nx*ny*k)]-bb[a+vars*(i+nx*(ny-3)+nx*ny*k)];
          }
          else{
            //Periodic
            xx[a+vars*(i+nx*ny*k)]=xx[a+vars*(i+nx*(ny-2)+nx*ny*k)];
            xx[a+vars*(i+nx*(ny-1)+nx*ny*k)]=xx[a+vars*(i+nx+nx*ny*k)];
            bb[a+vars*(i+nx*ny*k)]=bb[a+vars*(i+nx*(ny-2)+nx*ny*k)];
            bb[a+vars*(i+nx*(ny-1)+nx*ny*k)]=bb[a+vars*(i+nx+nx*ny*k)];
          }
        }
      }
    }

    //zBoundaries
    #pragma omp for
    for (int j = 1; j < ny-1; j++) {
      for (int i = 1; i < nx-1; i++) {
        for (int a = 0; a < vars; a++) {
          //z Boundary
          if (zbcon==1) {
            //Relfective
            xx[a+vars*(i+nx*j)]=xx[a+vars*(i+nx*j+nx*ny*2)];
            xx[a+vars*(i+nx*j+nx*ny*(nz-1))]=xx[a+vars*(i+nx*j+nx*ny*(nz-3))];
            bb[a+vars*(i+nx*j)]=bb[a+vars*(i+nx*j+nx*ny*2)];
            bb[a+vars*(i+nx*j+nx*ny*(nz-1))]=bb[a+vars*(i+nx*j+nx*ny*(nz-3))];
          }
          else if (zbcon==3) {
            //Outflow
            xx[a+vars*(i+nx*j)]=2*xx[a+vars*(i+nx*j+nx*ny)]-xx[a+vars*(i+nx*j+nx*ny*2)];
            xx[a+vars*(i+nx*j+nx*ny*(nz-1))]=2*xx[a+vars*(i+nx*j+nx*ny*(nz-2))]-xx[a+vars*(i+nx*j+nx*ny*(nz-3))];
            bb[a+vars*(i+nx*j)]=2*bb[a+vars*(i+nx*j+nx*ny)]-xx[a+vars*(i+nx*j+nx*ny*2)];
            bb[a+vars*(i+nx*j+nx*ny*(nz-1))]=2*bb[a+vars*(i+nx*j+nx*ny*(nz-2))]-bb[a+vars*(i+nx*j+nx*ny*(nz-3))];
          }
          else{
            //Periodic
            xx[a+vars*(i+nx*j)]=xx[a+vars*(i+nx*j+nx*ny*(nz-2))];
            xx[a+vars*(i+nx*j+nx*ny*(nz-1))]=xx[a+vars*(i+nx*j+nx*ny)];
            bb[a+vars*(i+nx*j)]=bb[a+vars*(i+nx*j+nx*ny*(nz-2))];
            bb[a+vars*(i+nx*j+nx*ny*(nz-1))]=bb[a+vars*(i+nx*j+nx*ny)];
          }
        }
      }
    }



    //###################################################
    //
    //      SECTION 4.2: Function computation.
    //      Each consituent equation is computed below
    //      the last, that being the element indicies of
    //      the equation must be N more than the last.
    //      The loops themselves are over the LOCAL indicies
    //      for that processor, so although the processor
    //      ownership range may begin midway through the
    //      vector in the global index this loop still needs
    //      to begin at 0.
    //
    //###################################################
    #pragma omp for
    for (int k=1; k<nz-1; k++) {
        for (int j=1; j<ny-1; j++) {
            for (int i=1; i<nx-1; i++) {

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



                //###################################################
                //
                //   SECTION 4.2.1: Calculate chi, the hall parameter
                //   set the magnetic field direction bdir and the
                //   transport coefficients
                //
                //####################################################

                //Calculate chi
                chi=hall(Z,bb,l); //this function calculates the hall parameter at mesh position i j k using the value of the magnetic field, density and temperature at that point

                //set magnetid field direction
                setbdir(l,bb,bdir);

                //set transport coefficients
                setcoeff(Z,chi,resistivity,bdir,1);
                setcoeff(Z,chi,conductivity,bdir,3);
                setcoeff(Z,chi,thermoelectric,bdir,2);

                //laser intensity values
                lasnorm=norm(laser[l/vars])+norm(laser2[l/vars])+2.0*theta*real(laser[l/vars]*conj(laser2[l/vars]));
                lasnormx1=norm(laser[(l+d1x)/vars])+norm(laser2[(l+d1x)/vars])+2.0*theta*real(laser[(l+d1x)/vars]*conj(laser2[(l+d1x)/vars]));
                lasnormx2=norm(laser[(l+d2x)/vars])+norm(laser2[(l+d2x)/vars])+2.0*theta*real(laser[(l+d2x)/vars]*conj(laser2[(l+d2x)/vars]));
                lasnormy1=norm(laser[(l+d1y)/vars])+norm(laser2[(l+d1y)/vars])+2.0*theta*real(laser[(l+d1y)/vars]*conj(laser2[(l+d1y)/vars]));
                lasnormy2=norm(laser[(l+d2y)/vars])+norm(laser2[(l+d2y)/vars])+2.0*theta*real(laser[(l+d2y)/vars]*conj(laser2[(l+d2y)/vars]));
                lasnormz1=norm(laser[(l+d1z)/vars])+norm(laser2[(l+d1z)/vars])+2.0*theta*real(laser[(l+d1z)/vars]*conj(laser2[(l+d1z)/vars]));
                lasnormz2=norm(laser[(l+d2z)/vars])+norm(laser2[(l+d2z)/vars])+2.0*theta*real(laser[(l+d2z)/vars]*conj(laser2[(l+d2z)/vars]));



                //###################################################
                //
                //      SECTION 4.2.3: Mesh interior, the function
                //      vector for each equation is now specified
                //      for the interior part of the mesh
                //
                //###################################################

                //Continuity equation
                ff[nin]=xx[nin]*(1.0+alpha*(xx[Vxin+d1x]-xx[Vxin+d2x])+beta*(xx[Vyin+d1y]-xx[Vyin+d2y])+gammad*(xx[Vzin+d1z]-xx[Vzin+d2z]))+alpha*xx[Vxin]*(xx[nin+d1x]-xx[nin+d2x])+beta*xx[Vyin]*(xx[nin+d1y]-xx[nin+d2y])+gammad*xx[Vzin]*(xx[nin+d1z]-xx[nin+d2z])-bb[nin];

                //X Momentum
                ff[Vxin]=R*xx[nin]*(xx[Vxin]-bb[Vxin]+alpha*xx[Vxin]*(xx[Vxin+d1x]-xx[Vxin+d2x])+beta*xx[Vyin]*(xx[Vxin+d1y]-xx[Vxin+d2y])+gammad*xx[Vzin]*(xx[Vxin+d1z]-xx[Vxin+d2z]))+0.5*(alpha*xx[nin]*(xx[Tin+d1x]-xx[Tin+d2x])+alpha*xx[Tin]*(xx[nin+d1x]-xx[nin+d2x]))-dt*xx[jyin]*xx[Bzin]+dt*xx[jzin]*xx[Byin]+alpha*bb[nin]*(lasnormx1-lasnormx2)*Pc;

                //Y Momentum
                ff[Vyin]=R*xx[nin]*(xx[Vyin]-bb[Vyin]+alpha*xx[Vxin]*(xx[Vyin+d1x]-xx[Vyin+d2x])+beta*xx[Vyin]*(xx[Vyin+d1y]-xx[Vyin+d2y])+gammad*xx[Vzin]*(xx[Vyin+d1z]-xx[Vyin+d2z]))+0.5*(beta*xx[nin]*(xx[Tin+d1y]-xx[Tin+d2y])+beta*xx[Tin]*(xx[nin+d1y]-xx[nin+d2y]))+dt*xx[jxin]*xx[Bzin]-dt*xx[jzin]*xx[Bxin]+beta*bb[nin]*(lasnormy1-lasnormy2)*Pc;

                //Z Momentum
                ff[Vzin]=R*xx[nin]*(xx[Vzin]-bb[Vzin]+alpha*xx[Vxin]*(xx[Vzin+d1x]-xx[Vzin+d2x])+beta*xx[Vyin]*(xx[Vzin+d1y]-xx[Vzin+d2y])+gammad*xx[Vzin]*(xx[Vzin+d1z]-xx[Vzin+d2z]))+0.5*(gammad*xx[nin]*(xx[Tin+d1z]-xx[Tin+d2z])+gammad*xx[Tin]*(xx[nin+d1z]-xx[nin+d2z]))-dt*xx[jxin]*xx[Byin]+dt*xx[jyin]*xx[Bxin]+gammad*bb[nin]*(lasnormz1-lasnormz2)*Pc;

                //Energy
                ff[Tin]=1.5*xx[nin]*(xx[Tin]-bb[Tin]+alpha*xx[Vxin]*(xx[Tin+d1x]-xx[Tin+d2x])+beta*xx[Vyin]*(xx[Tin+d1y]-xx[Tin+d2y])+gammad*xx[Vzin]*(xx[Tin+d1z]-xx[Tin+d2z]))+2*((alpha*(xx[qxin+d1x]-xx[qxin+d2x])+beta*(xx[qyin+d1y]-xx[qyin+d2y])+gammad*(xx[qzin+d1z]-xx[qzin+d2z])))+xx[nin]*xx[Tin]*(alpha*(xx[Vxin+d1x]-xx[Vxin+d2x])+beta*(xx[Vyin+d1y]-xx[Vyin+d2y])+gammad*(xx[Vzin+d1z]-xx[Vzin+d2z]))-2*dt*((xx[Exin]+xx[Vyin]*xx[Bzin]-xx[Vzin]*xx[Byin])*xx[jxin]+(xx[Eyin]+xx[Vzin]*xx[Bxin]-xx[Vxin]*xx[Bzin])*xx[jyin]+(xx[Ezin]+xx[Vxin]*xx[Byin]-xx[Vyin]*xx[Bxin])*xx[jzin])-dt*IBc*bb[nin]*bb[nin]*(lasnorm/pow(bb[Tin],1.5));

                //X Ohms law
                bierx=0;//0.5*alpha*(xx[Tin]*(xx[nin+d1x]-xx[nin+d2x])+xx[nin]*(xx[Tin+d1x]-xx[Tin+d2x]))/dt;
                ff[Exin]=xx[nin]*(xx[Exin]+xx[Vyin]*xx[Bzin]-xx[Vzin]*xx[Byin])+bierx-(resistivity[0]*xx[jxin]+resistivity[1]*xx[jyin]+resistivity[2]*xx[jzin])/pow(xx[Tin],1.5)+0.5*xx[nin]*(thermoelectric[0]*alpha*(xx[Tin+d1x]-xx[Tin+d2x])+thermoelectric[1]*beta*(xx[Tin+d1y]-xx[Tin+d2y])+thermoelectric[2]*gammad*(xx[Tin+d1z]-xx[Tin+d2z]));

                //Y Ohms law
                biery=0;//0.5*beta*(xx[Tin]*(xx[nin+d1y]-xx[nin+d2y])+xx[nin]*(xx[Tin+d1y]-xx[Tin+d2y]))/dt;
                ff[Eyin]=xx[nin]*(xx[Eyin]+xx[Vzin]*xx[Bxin]-xx[Vxin]*xx[Bzin])+biery-(resistivity[3]*xx[jxin]+resistivity[4]*xx[jyin]+resistivity[5]*xx[jzin])/pow(xx[Tin],1.5)+0.5*xx[nin]*(thermoelectric[3]*alpha*(xx[Tin+d1x]-xx[Tin+d2x])+thermoelectric[4]*beta*(xx[Tin+d1y]-xx[Tin+d2y])+thermoelectric[5]*gammad*(xx[Tin+d1z]-xx[Tin+d2z]));

                //Z Ohms law
                bierz=0;//0.5*gammad*(xx[Tin]*(xx[nin+d1z]-xx[nin+d2z])+xx[nin]*(xx[Tin+d1z]-xx[Tin+d2z]))/dt;
                ff[Ezin]=xx[nin]*(xx[Ezin]+xx[Vxin]*xx[Byin]-xx[Vyin]*xx[Bxin])+bierz-(resistivity[6]*xx[jxin]+resistivity[7]*xx[jyin]+resistivity[8]*xx[jzin])/pow(xx[Tin],1.5)+0.5*xx[nin]*(thermoelectric[6]*alpha*(xx[Tin+d1x]-xx[Tin+d2x])+thermoelectric[7]*beta*(xx[Tin+d1y]-xx[Tin+d2y])+thermoelectric[8]*gammad*(xx[Tin+d1z]-xx[Tin+d2z]));

                //X Faraday
                ff[Bxin]=xx[Bxin]-bb[Bxin]+beta*(xx[Ezin+d1y]-xx[Ezin+d2y])-gammad*(xx[Eyin+d1z]-xx[Eyin+d2z]);

                //Y Faraday
                ff[Byin]=xx[Byin]-bb[Byin]+gammad*(xx[Exin+d1z]-xx[Exin+d2z])-alpha*(xx[Ezin+d1x]-xx[Ezin+d2x]);

                //Z Faraday
                ff[Bzin]=xx[Bzin]-bb[Bzin]+alpha*(xx[Eyin+d1x]-xx[Eyin+d2x])-beta*(xx[Exin+d1y]-xx[Exin+d2y]);

                //X Ampere
                ff[jxin]=dt*mu*xx[jxin]+gammad*(xx[Byin+d1z]-xx[Byin+d2z])-beta*(xx[Bzin+d1y]-xx[Bzin+d2y]);

                //Y Ampere
                ff[jyin]=dt*mu*xx[jyin]+alpha*(xx[Bzin+d1x]-xx[Bzin+d2x])-gammad*(xx[Bxin+d1z]-xx[Bxin+d2z]);

                //Z Ampere
                ff[jzin]=dt*mu*xx[jzin]+beta*(xx[Bxin+d1y]-xx[Bxin+d2y])-alpha*(xx[Byin+d1x]-xx[Byin+d2x]);

                //limiting procedure
                qlimx=0.25*pow(xx[Tin],2.5)*(conductivity[0]*alpha*(xx[Tin+d1x]-xx[Tin+d2x])+conductivity[1]*beta*(xx[Tin+d1y]-xx[Tin+d2y])+conductivity[2]*gammad*(xx[Tin+d1z]-xx[Tin+d2z]))/dt+0.5*xx[Tin]*(thermoelectric[0]*xx[jxin]+thermoelectric[1]*xx[jyin]+thermoelectric[2]*xx[jzin]);
                qlimy=0.25*pow(xx[Tin],2.5)*(conductivity[3]*alpha*(xx[Tin+d1x]-xx[Tin+d2x])+conductivity[4]*beta*(xx[Tin+d1y]-xx[Tin+d2y])+conductivity[5]*gammad*(xx[Tin+d1z]-xx[Tin+d2z]))/dt+0.5*xx[Tin]*(thermoelectric[3]*xx[jxin]+thermoelectric[4]*xx[jyin]+thermoelectric[5]*xx[jzin]);
                qlimz=0.25*pow(xx[Tin],2.5)*(conductivity[6]*alpha*(xx[Tin+d1x]-xx[Tin+d2x])+conductivity[7]*beta*(xx[Tin+d1y]-xx[Tin+d2y])+conductivity[8]*gammad*(xx[Tin+d1z]-xx[Tin+d2z]))/dt+0.5*xx[Tin]*(thermoelectric[6]*xx[jxin]+thermoelectric[7]*xx[jyin]+thermoelectric[8]*xx[jzin]);

                qmag=sqrt(pow(qlimx,2.0)+pow(qlimy,2.0)+pow(qlimz,2.0));
                qfs=flim*bb[nin]*bb[Tin];
                  //X Heat flow
                ff[qxin]=xx[qxin]+qlimx/(1.0+qmag/qfs);

                  //Y Heat flow
                ff[qyin]=xx[qyin]+qlimy/(1.0+qmag/qfs);

                  //Z Heat flow
                ff[qzin]=xx[qzin]+qlimz/(1.0+qmag/qfs);

            }
        }
    }
    }



    /* Restore vectors */
    VecRestoreArray(ref,&xx);
    VecRestoreArray(b,&bb);
    VecRestoreArray(f,&ff);
    VecDestroy(&ref);
    return 0;
}

//##############################################################################
//
//      SECTION 5: Form the initial guess of the vector x, this must be done
//      or zeros will turn up in the diagonals of the jacobian leading to a
//      failed solver
//
//##############################################################################


PetscErrorCode FormInitialGuess(SNES snes, Vec x,void *ctx){

    datastruct       *user = (datastruct*) ctx;

    VecAYPX(x,0,user->b);
    return 0;
}

//##############################################################################
//
//      SECTION 6: Paraxial laser solver
//      input: data structure ,
//
//##############################################################################

PetscErrorCode PARAXIAL(datastruct &ctx, complex<double> *laser){

  Vec           b=ctx.b;
  //complex<double> *laser=ctx.laser;
  PetscScalar  *bb;

  complex<double> diagx[nx];
  complex<double> bx[nx];
  complex<double> diagy[ny];
  complex<double> by[ny];

  VecGetArray(b,&bb);

  complex<double> A1,A2,A3,A4;
  int q;
  double ep0=1.0-1/nc; //dielectric constant
  double re1,re2, im1,im2,dens,temp;


  //#################################
  //
  //    Solve the laser field
  //
  //#################################

  for (int k = 1; k < nz-1; k++) {


    //yboundary
    for (int i = 0; i < nx; i++) {
      laser[i+nx*ny*k]=laser[i+nx*(ny-2)+nx*ny*k];
      laser[i+nx*(ny-1)+nx*ny*k]=laser[i+nx+nx*ny*k];
    }



    A2=complex<double>(dz/dx/dx,0.0);
    A4=complex<double>(-dz/dy/dy,0.0);

    //First Sweep in y direction
    for (int j = 1; j < ny-1; j++) {

      //form vectors for linear solve in x
      for (int i = 1; i < nx-1; i++) {

        q=(i+nx*j+nx*ny*k);

        dens=bb[q*vars];
        temp=bb[q*vars+4];
        re1=-2.0*dz/dx/dx+wn*wn*dz*0.5*(1-dens)/nc;
        im1=4.0*wn+wn*wn*dz*0.5*dens*dens/(nc*omega*tei0*pow(temp,1.5));
        re2=2.0*dz/dy/dy-wn*wn*dz*0.5*(1-dens)/nc;
        im2=4.0*wn-wn*wn*dz*0.5*dens*dens/(nc*omega*tei0*pow(temp,1.5));

        A1=complex<double>(re1,im1);
        A3=complex<double>(re2,im2);

        //diagonal
        diagx[i-1]=A1;
        //rhs
        bx[i-1]=A4*(laser[q+nx]+laser[q-nx])+A3*laser[q];

      }
      //linear solve
      SM(A2,diagx,A2,bx,nx-2);

      //Reset values in u
      for (int i = 1; i < nx-1; i++) {
        q=(i+nx*j+nx*ny*(k+1));
        laser[q]=bx[i-1];
      }

    }


    //xboundary
    for (int j = 0; j < ny; j++) {
      laser[nx*j+nx*ny*(k+1)]=laser[nx-2+nx*j+nx*ny*(k+1)];
      laser[nx-1+nx*j+nx*ny*(k+1)]=laser[1+nx*j+nx*ny*(k+1)];
    }


    A2=complex<double>(dz/dy/dy,0.0);
    A4=complex<double>(-dz/dx/dx,0.0);
    //Sweep in x direction

    for (int i = 1; i < nx-1; i++) {

      //form vectors for linear solve in y
      for (int j = 1; j < ny-1; j++) {

        q=(i+nx*j+nx*ny*(k+1));

        dens=bb[vars*(q-nx*ny)];
        temp=bb[vars*(q-nx*ny)+4];
        re1=-2.0*dz/dy/dy+wn*wn*dz*0.5*(1-dens)/nc;
        im1=4.0*wn+wn*wn*dz*0.5*dens*dens/(nc*omega*tei0*pow(temp,1.5));
        re2=2.0*dz/dx/dx-wn*wn*dz*0.5*(1-dens)/nc;
        im2=4.0*wn-wn*wn*dz*0.5*dens*dens/(nc*omega*tei0*pow(temp,1.5));

        A1=complex<double>(re1,im1);
        A3=complex<double>(re2,im2);

        //diagonal
        diagy[j-1]=A1;
        //rhs
        by[j-1]=A4*(laser[q+1]+laser[q-1])+A3*laser[q];

      }
      //linear solve
      SM(A2,diagy,A2,by,ny-2);

      //Reset values in u
      for (int j = 1; j < ny-1; j++) {
        q=(i+nx*j+nx*ny*(k+1));
        laser[q]=by[j-1];
      }

    }

  }


  VecRestoreArray(b,&bb);

  return 0;
}

//##############################################################################
//
//      SECTION 7: complex tridiagonal solver
//      input: subdiagonal term, diagonal vector, super diagonal, RHS vector and length
//
//##############################################################################

void zTriSolver(complex<double> subdiag, complex<double> *diag, complex<double> superdiag, complex<double> *b,int length){

  //Tridiagonal matrix solver with Sherman-Morrison correction for periodic boundary terms
  complex<double> arr1[length];
  complex<double> arr2[length];

  //Forward
  arr1[0]=diag[0];
  arr2[0]=b[0]/diag[0];

  for (int i = 1; i < length; i++) {
    arr1[i]=diag[i]-subdiag*superdiag/arr1[i-1];
    arr2[i]=(b[i]-subdiag*arr2[i-1])/arr1[i];
  }
  //Backward
  b[length-1]=arr2[length-1];
  for (int i = length-2; i >= 0; i--) {
    b[i]=arr2[i]-superdiag*b[i+1]/arr1[i];
  }


}

//##############################################################################
//
//      SECTION 8: Sherman morrison solver for periodic/cyclic tridiagonal
//      input: subdiagonal term, diagonal vector, super diagonal, RHS vector and length
//
//##############################################################################


void SM(complex<double> subdiag, complex<double> *diag, complex<double> superdiag, complex<double> *b,int length){

  complex<double> gamma,lambda;
  gamma=b[0];
  diag[0]+=gamma;
  diag[length-1]+=superdiag*subdiag/gamma;

  zTriSolver(subdiag,diag,superdiag,b,length);

  complex<double> u[length];
  for (int i = 0; i < length; i++) {
    u[i]=complex<double>(0.0,0.0);
  }
  u[0]=-gamma;
  u[length-1]=subdiag;

  zTriSolver(subdiag,diag,superdiag,u,length);

  lambda=(b[0]-subdiag*b[length-1]/gamma)/(1.0+u[0]-subdiag*u[length-1]/gamma);
  for (int i = 0; i < length; i++) {
    b[i]-=lambda*u[i];
  }


}
