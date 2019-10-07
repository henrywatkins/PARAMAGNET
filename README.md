## PARAMAGNET
Long-Pulse Laser-Plasma physics simulation code
-----------------------------------------------

This code has been written to simulate laser-plasma interactions in the long (~1ns) pulse regime. Given the collisionality of this regime allows us to model the plasma using MagnetoHydrodynamics (MHD) coupled to a Paraxial laser solver.   

Prerequisites:
--------------

- [PETSc](https://www.mcs.anl.gov/petsc/) - this has been tested with version 3.8
- [OpenMP](https://www.openmp.org/)

How to install/compile:
-----------------------

Clone the repository and cd into the main directory. Ensure the environment variables are set for the PETSc directory; this should be in the .bashrc file. Install with the command `make all`. This will create a binary in the bin directory.

How to define a simulation problem:
-----------------------------------

Modify the files Initial.cpp and Global.cpp to set the initial condition and simulation parameters respectively. 

The Global.cpp file contains global variables that are set at compile time. 
The `Physical State` section contains definitions of:
* the reference density `n0`,
* the reference temperature `T0`, 
* the ionization state `Z`(which at present only works up to Z=8),
* the period of the simulation (in seconds)
* the physical dimensions of the simulation box (x,y,z, in metres)
The `Collisional Parameters` should not be changed.
The `Mesh Variables` section contains a set of parameters that set the mesh parameters. The only ones that should be changed are the variables `m,nx,ny,nz` which set the number of timesteps and the number of cells in x,y,z respectively.
The `Laser Variables` section contains the parameters for the laser solver. The parameters one can change are:
* The wavelength `wl` in metres,
* The intensity `Intensity` in watts per square metre,
The `Output Varibales` section defines whether a field should be saved to file as output from the simulation, and `tfrac` defines how many steps there are between output dumps (`tfrac=10` means output is dumped every 10 timesteps).
The `Boundary conditions` section defines the boundaries in x, y and z.

In the Inital.cpp file one can set the initial profile of the various fields. To set a profile, find the discrete form of the profile on a cartesian mesh (with cell position i,j,k for the x,y and z dimensions respectively). Ensure the profile is set in normalised units. The `initvalue` array contains the intial state of the code, with the indicies `nin` representing number density, `Tin` temperature, `Vxin` the x velocity and `Bxin` the magnetic field in x. The laser field is defined on the `z=0` boundary, as such only profiles with `k=0` should be used. 

Once all modifications to the Global and Initial files are complete, recompile the source using `make all`.

How to run:
-----------

Set the number of OpenMP threads using `OPENMP_NUM_THREADS=N` where the number of threads is N. Run with the command `./bin/MAGNET`. The output can be found in raw ASCII format and VTK format in the `output` directory.
