Output
======

The PARAMAGNET code produces two kinds of output, plain ASCII text files
for each variable that is output by the code and VTK output files.

ASCII Output
------------
The plain ASCII text file output can be used for analysis in Matlab or Python.
There is a file for each variable, density, temperature etc.

VTK Output
----------
The VTK file output can be used for analysis in Paraview or Visit.
This is useful for 3D visualisation and to quickly look at the output of the code.
The format used for the VTK output is the legacy .vtk ASCII format.
It works and was simple to code up but more modern and sophisticated
formats exist in the VTK family.
