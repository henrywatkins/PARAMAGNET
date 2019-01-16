Simulation Setup
================

Initialisation
--------------
To setup a simulation changes must be made to the ``Initial.cpp`` file and
the ``Global.cpp`` file. The ``Global.cpp`` file contains the globally defined
variables in the code.

The ``Global.cpp`` file contains:

  * Physical State: variables where the temperature, density, domain size, time period and ionisation state is set
  * Collisional Parameters: these should not be changed
  * Mesh Variables: change the number of timesteps and mesh size
  * Laser Variables: change the peak intensity and wavelength, also one can turn off the ponderomotive and inverse bremmsstrahlung terms if needed
  * Output Variables: change the variables to output and how many timesteps should be Output
  * Boundary Conditions: choose between periodic and outflow boundary conditions

The ``Initial.cpp`` file contains:

  * Initial Condition: this sets the initial state of the plasma in normalised units.
  * Initial Laser: this sets the laser profile at the z=0 boundary

OpenMP
------
To run in parallel one first has to specify the number of OpenMP threads
that the program will use. Set this with the command:
``export OMP_NUM_THREADS=N``
with N as the number of shared memory processes needed.

Compilation
-----------
To compile the code run the command:
``make MAGNET``
The makefile should take care of the compilation and linking process.
Within the makefile one must set the compiler variable to a C++ compiler
available on the computer.

Running
-------
To run the code use the command:
``./MAGNET``
