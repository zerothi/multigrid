## Multigrid method for Poisson equation

This simple code solves the Poisson equation for open boundary problems
often encountered in DFT+NEGF calculations.

In particular can it solve problems with periodic, Dirichlet and Neumann
boundary conditions individually selected on each boundary.

### Method

We employ an adaptive/customizable multi grid method where every level
can be controlled arbitrarily by the user.

### Usage

Compile the MG library using

    make

to compile the executable for solving input files, do

	make mg

which creates the executable `src/mg` which can read input files.
In the `test` folder several examples exists of input files for
complex MG setups.

#### Precision

The precision can be controlled to arbitrary value. However, it is not
stringently following conventional multi-grid methods as we simply want
a reasonable guess for the Poisson solution.
