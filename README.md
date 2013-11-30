IBM-ecology-simulator
=====================

Individual Base Model simulator. C library that allow to run eficient IBM simulations with Monte Carlo algorithms

PDSpaceStructure.c is the driver or main file where the main function is.

libPP_5.0.c is the source code of all the important functions. It is system indepent.

EntSalArb_MP.c the source code that write the structures defined in libPP_X.c to hard drive. Just a function to 
create a directory is system dependent (linux).

GNA.c is the random number generator used by libPP_X.c



