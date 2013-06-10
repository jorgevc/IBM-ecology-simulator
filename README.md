IBM-ecology-simulator
=====================

Individual Base Model simulator. C library that allow to run eficient IBM simulations with Monte Carlo algorithms

PDSpaceStructure.c is the driver or main file where the main function is.

libPP_5.0.c is the source code of all the important functions. It is system indepent.

EntSalArb_MP.c the source code that write the structures defined in libPP_X.c to hard drive. Just a function to 
create a directory is system dependent (linux).

GNA.c is the random number generator used by libPP_X.c

You are not going to is this (at least at first) I don't plan to document it (well, just if I have some petitions):
ControlDinamico.c is a template to implement a real time interface for the simulations.

CD.c This just work if ControlDinamico.c is properly set in PDSpaceStructure.c . The client that comunicates in real time with the simulation: stop the process, write results to disc, 
show step iteration, pause, set a crono to stop, write a distribution to disc.  I implement this to see the results of a 
running simulation without stopping it and decide if it should continue or it can be stopped. 

