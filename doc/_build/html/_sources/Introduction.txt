Introduction
============

The PARAMAGNET is a 3D plasma physics code for the simulation of magnetized
collisional plasmas with full electron transport in a Cartesian geometry.

Physics
-------
The code uses a MHD model for the plasma with the full Braginskii electron
transport coefficients for the anisotropic heat flow, resitivity and
thermoelectric coefficient.

This plasma model is coupled to a paraxial laser model via inverse
bremmsstrahlung and the ponderomotive force.

These two modules are separate and can be solved independently.

Numerical Methods
-----------------
The plasma solver uses a fully implicit finite difference scheme.
This produces a large nonlinear matrix problem that is solved using
the Jacobian-Free Newton method with a GMRES linear solver provided
by the PETSc library.

The laser solver uses finite differences. This produces a
sweeping across the z direction where a linear problem has
to be solved at each slice of the domain along z. To solve this
linear problem the Alternating-Direction-Implicit method is used in
compination with a direct complex tridiagonal matrix solver.
