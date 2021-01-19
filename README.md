# 2D_Diffusion
> Bidimensional Solver of Diffusion Equation with Finite Differences - FORTRAN

2D_Difussion finite differences solver with CSC storage and solvers for matrices (Thematic) or ADI solver(Developing State)

With Fortran - GNU-Fortran compiler (and Matlab originals).

## Installation

OS X & GNU-Linux:

```sh
git clone https://github.com/sacastiblancob/2D_Diffusion.git
```

```sh
make clean; make
```

## Usage example

To run the solver you should have a Fotran compiler (if possible use GNU-Fortran!!), change in the Makefile the lines related with the compiler and user specific configuration (the actual version works perfectly with Gfortran).

Once the Makefile is modified (if you need to), you should compile it by typing the make clean and make commands. The .out executable file will be located in the ./bin/ folder.

Modify the configuration in the ./diffconf.nml file, if you want to change some default variable, modifiable variables are explained within this file.

Once you have configured your scenario and compile the solver, you can run it by the following command:

```sh
./bin/diff_df.out
```

Wait to reach the Final Time and all the information you need to know will be printed in the terminal, it will show you values such as the tolerance and number of iterations used by the Conjugate Gradient Solver with CSC storage for matrices.

The results will be storaged in the ./res/ file, within the .dat files. You can load them into Paraview by using the "TecPlot Reader" option.

Printed Variables: Tracer (C), Solution Abs Error (ERR)

## Folder Contents

./bin/ --> Executable file

./matlab/ --> Matlab Files (main = Solver_PF.m)

./mod/ --> Fortran modules

./obj/ --> Fortran objects

./res/ --> Solver results(.dat files)

./src/ --> Fotran source codes

./diffconf.nml --> Configuration file

./Makefile --> Make compiler file

./README.md --> You are standing here

## Meta
Sergio A. Castiblanco-Ballesteros
Bogota - Colombia
sergio.castiblanco@javeriana.edu.co
sacastiblancob@unal.edu.co

> Free Distribution and Open Source (As it should be!!)


