# TWI2D-Code

Date: February 20th, 2023.
Author: Júlio Oliva Frigério
Author's email: juliojof@sep.stanford.edu

License: GNU Lesser General Public License v2.1


First TWI2D code repository. In this README file you will find esssential instructions to compile and run the software. But to execute the software you should also rely on the Jupyter Notebooks provided in the sister repository "TWI2D-JupyterNotebooks".


=========================
====== Compilation ======
=========================

** Requirements **
To compile this code it is necessary to have the following software:
- Cuda. Any version that is not too early will do, but you have to change the Makefile so that it links with the correct version.
- FFTW3.
- MPI with compiler mpicc.
- Compiler g++.

** Instructions **
To compile this code you just need to run make from within the directory where the code (.c and .h files) is located. The compilation should conclude without errors. After that the code is compiled.



=========================
======= Execution =======
=========================

** Parameters file **
The software takes as input a file with execution parameters. This file must be named InputParms and must be located in the folder where the software 
is being executed. Templates of this file are provided in this repository. It is also possible to build these parameters files with the Jupyter Notebooks provided in the sister repository "TWI2D-JupyterNotebooks".

** Running instructions **
To execute (run) the software you should run the following command:
> mpirun -n ## PATH/TO/EXECUTABLE
In the command above, ## is an integer number that defines the number of processes executed by the MPI.
