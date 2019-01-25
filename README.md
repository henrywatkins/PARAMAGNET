# PARAMAGNET
Long-Pulse Laser-Plasma physics simulation code
-----------------------------------------------
Prerequisites:

	- [PETSc](https://www.mcs.anl.gov/petsc/) - this has been tested with version 3.7
	- [OpenMP](https://www.openmp.org/)

How to install/compile:
Clone the repository and cd into the main directory. Ensure the environment variables are set for the PETSc directory; this should be in the .bashrc file. Install with the command

	`make all`

This will create a binary in the bin directory.

How to define a simulation problem:
Modify the files Initial.cpp and Global.cpp to set the initial condition and simulation parameters respectively. Recompile the source using `make all` as above.

How to run:
Set the number of OpenMP threads using `OPENMP_NUM_THREADS=N` where the number of threads is N. Run with the command `./bin/MAGNET`. The output can be found in raw ASCII format and VTK format in the output directory. 
