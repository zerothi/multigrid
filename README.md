## Multigrid method for Poisson equation

This simple code solves the Poisson equation for open boundary problems
often encountered in DFT+NEGF calculations.

In particular can it solve problems with periodic, Dirichlet and Neumann
boundary conditions individually selected on each boundary.

### Method

We employ an adaptive/customizable multi grid method where every level
can be controlled arbitrarily by the user.

#### Precision

The precision can be controlled to arbitrary value. However, it is not
stringently following conventional multi-grid methods as we simply want
a reasonable guess for the Poisson solution.
