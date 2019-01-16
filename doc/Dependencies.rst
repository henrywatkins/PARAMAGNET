Dependencies
============

PETSc
-----
PETSc is a library of linear and nonlinear solvers that are required to run
the PARAMAGET code. The core solvers lie at the heart of the plasma solver
within PARAMAGNET and so must be installed. In order to do so go to the PETSc
website https://www.mcs.anl.gov/petsc/ and follow the download and installation
instructions. The version used in the code was **3.6.3** and so it is best to
download this version.

OpenMP
------
The PARAMAGNET code uses OpenMP parallelisation for parts of the code and so a
version of OpenMP must also be installed. The download and installation
instructions can be found on the website https://http://www.openmp.org/.
OpenMP is a shared memory model parallel library and so it can only be used
for a single node. 
