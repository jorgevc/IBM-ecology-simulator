------------
Dependecies 
----------
-Install fftw3 library from http://www.fftw.org/download.html
-openMP http://openmp.org in gcc it is already ready for use.
-gnuplot is needed if quick visualitation is needed http://sourceforge.net/projects/gnuplot/

-----------
Compilation
-----------
Once fftwe is installed the
Program can be built using default make arguments.
i.e. just type "make" in the command line

---------
Run the simulations
---------
Once compiled
run ./crece.out to run the actual simulation
The following files will be created:

a) A set of size distributions at diferent time steps
	./DATOS_TAM/DT_#	where # stands for a time step of the simulation (units of time ~ 10 weeks)
b) A set of field representations at diferent time steps (units of time ~ 10 weeks)
	./DATOS_TAM/T_#		where # is the time step of the simulation (units of time ~ 10 weeks)
c) A table containing the instant growth rate
	./DATOS_TAM/GrowthR
d) A table containing the instant resource intake
	./DATOS_TAM/ResourceR
e) A table containig the mean size, mean density vs time
	./DATOS_TAM/thinning
	
---------
Some visualitations
---------
Once the programm is finished
run:
gnuplot "plots.gp"   (this can take some minutes)
in order to generate the following visualitations
	./dist.gif	animation of the size distribution (ensemble average)
	./field.gif	animation of one simulated field
	./MeanDensity.pdf	Mean density vs time (ensemble average)
	./MeanSize.pdf		Mean size vs time (ensemble average)
	./ResourceConsumption	Mean resource consumption vs time (ensemble average)

------------------------------
Important Files and functions
------------------------------
-Competencia.c:
This is the driver of the program where the main function is.
The general structure of the main function is
1) Initialization of the variables
2) Memory allocation for the tracking sistem (arrays where diferent rates are going to be stored)
3) Parallel code using OMP:
	3.1) Memory allocation for the ensemble of simulations
	3.2) Setting the initial spatial distribution of individuals
	3.3) Initialization of thread traking sistem (to improve performace by avoiding fake memory sharing)
	3.4) Iteration of each time steep. (Each time steep performs a Monte Carlo sweep)
		3.4.1) Loop through the simulations that correspond to "this" thread
		3.4.2) Keep track of mean density, mean size and map between MC sweep and fisical time
		3.4.3) After certain number of steeps write to disk by coping the distribution of each simualtion
		 to a shared memory array and averaging. Then write to disk
	3.5) Copy tracking variables of each thread to a shared memory (global) tracking sistem and then Free memory
4) Write tracking variables to disk and free remaining memory

-libPP_6.1.c, libPP_6.1.h:
This is the main library. It contains many helper functions but the core two are
1) BarrMCcRyCampTamano
	this function iterates unit events until a monte carlo sweep is performed.
	Here are mesured the instant rates
2)ActualizaUniv
	This is the unit event of the poisson process in the MC simulation
	Here is where the model is actually implemented
	
-EntSalArb_MP_Comp.c, EntSalArb_MP_Comp.h
This library contain helper functions to write the spetial structures to disk

-GNA.c , GNA.h
This is the implementation of the random number generator (KISS)

---------
Licence
---------
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */
