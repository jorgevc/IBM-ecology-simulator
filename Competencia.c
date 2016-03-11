/*
Copyright 2012 Jorge Velazquez
*/
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "libPP_6.1.h"
#include "EntSalArb_MP_Comp.h"
#include "GNA.h"
#include <math.h>


main(){	
#ifdef VIRTUAL_GRID			///< The algorithm should be define at top of file "libPP_6.1.h"
printf("VIRUAL_GRID=TRUE\n");		///< Prints the algorithm being used "Virtual grid of resources" deprecated
#endif
#ifdef SOI
printf("SOI=TRUE\n");		///< Prints the algorithm being used "Zone of Influence"
#endif

/** This section initialize the simulation parameteres */
runDescriptor run;		///< structure containing the simulation parameters
run.X=200; 			 	///< X lenght of the simulated field
run.Y=run.X;			///< Y lenght of the simulated field
run.grid_units=1.0; 	///< convertion factor between lenght units to virtual grid cell lenght.
run.size_units=1.0; 	///< convertion factor between lenght units and virtual tree size units.
run.T_max=3001;			///< Number of monte carlo sweeps in the simulation
run.NoEnsambles=700;	///< Number of simulations in the ensemble of simulations.
int const write_interval=100;	///< Number of monte carlo sweeps between writting state to disk.

/** Basic model parameters **/
float a = 0.0415;			///< convertion factor between m^{3/4} and resource use 
float b = 0.0013;			///< Resource cost for maintenance per unit biomass per time 
run.Model.competitionAsymetry = 10.0;	///< The asymetry of competition "p".
///************************///

/// in this case the Area_units and Length_units are units of the imaginary resources grid.
/// Some convertion factors between computer simualation units and real fisical units
run.Model.ResourcesScale = 1.0; ///< convertion factor between lenght units of resource grid and simulation length.
a = a*0.19;		///< convertion to units of time of 10 weeks (not years)
b = b*0.19;		///< convertion to units of time of 10 weeks
double const Area_units=(run.grid_units*run.Model.ResourcesScale)*(run.grid_units*run.Model.ResourcesScale);
double const Length_units=run.grid_units*run.Model.ResourcesScale;

double coagulation_units = 1.0;
#ifdef SOI
coagulation_units = 52.0;
run.Model.coagulation_factor = coagulation_units*1.0;
#else
run.Model.coagulation_factor = coagulation_units*1.0*3.1416;
#endif

///* Some convertion factors to run the simulation of the evolution of the canopy diameter directly *///
run.Model.coagulation_radio_exp=1.0;  ///< Alometric relation: exponent relating steem radio and canopy radio 
run.Model.coagulation_exp=2*run.Model.coagulation_radio_exp; 
run.Model.Cr= 4.7;		///< convertion factor between m^{3/8} and canopy radio units
run.Model.coagulation_radio_factor= 1.0;
run.Model.metabolic_exp=2.6666; ///< Alometric relation: exponent relating canopy diameter and biomass (8/3)
run.Model.health_factor=0.0; ///< Fraction of allowed missed biomass before individual dies.
run.Model.Cg= (3.1416*run.Model.Cr*run.Model.Cr)/(a); ///< Resources needed for unit size increase
run.Model.Cm= run.Model.Cg*b;    ///<  Resources needed for metabolic needs
run.Model.growth_constant=(3.0*pow(run.Model.Cr,8.0/3.0))/(pow(2.0,1.0/3.0)*run.Model.Cg); ///<  needed resources per unit size increse.
run.Model.resource_rate=1.0/coagulation_units; 
run.initialMeanDistance=run.grid_units*1;  ///< mean distance between individuals (when random pattern is generated)
run.initialMinSeparation=1;	///< Minimun separation between individual in a random pattern (filter)
run.Model.birth_rate=0.0;	///< birth rate of new individuals
run.Model.RadioBirth=15;	///< radius of new offsprings from parent 
run.Model.dead_rate=0.0;	///< instrinsic dead rate of individuals
run.Model.intra_coagulation=0.0; ///< deprecated
///* Initial values for the trees *///
Individual indv;
indv.species=1;			///< species type
indv.size_float=2.5;	///< biomass of new individuals floating point variable
indv.size=run.size_units*10.0*indv.size_float;
indv.radio_float=R(indv, (&run.Model));  ///< dbh radio of individuals
indv.radio=indv.radio_float;	
indv.metabolism=0;	///< deprecated, used to improve time step performance
indv.health=0; 		///< deprecated, used to track individual healt
omp_set_num_threads(4);		///< sets the number of threads used by omp

/** internal initialization variables **/
float const separation=(float)run.initialMeanDistance;
int T_max = run.T_max;
int NoEnsambles=run.NoEnsambles;
run.Model.meta_needs=SetMetaNeeds(run.Model, Area_units); 
run.Model.R=SetR(run.Model, Length_units);
run.Model.M=SetM(run.Model, Area_units);
model modelo;
modelo = run.Model;
int NDX=(run.X*run.grid_units)+1;
int NDY=(run.Y*run.grid_units)+1;
double delta_s=1.0/(run.size_units);
///************************************************///

/** Initialization of sturctures that keep the resoults **/
	
Float1D_MP TamDist_1;		///< Instance of Float1D_MP that will store the size distribution at each steep
		InicializaFloat1D_MP(&TamDist_1, T_max+10);		///< initialization function of Float1D structures.
		
Float1D_MP CumulativeTamDist_1;		///< instance of Float1D_MP that can store the comulative size distribution
		InicializaFloat1D_MP(&CumulativeTamDist_1, T_max+10);
		
Float1D_MP MP_CorrelacionG;			///< instance of Float1D_MP that can store a pair correlation 
	InicializaFloat1D_MP(&MP_CorrelacionG, NDX);  ///<  initialization
	ResetFloat1D_MP(&MP_CorrelacionG);		///< we ensure it is full of zeros

char contenedor[150];					
	sprintf(contenedor,"DATOS_TAM");		///< Name of the directory where data is going to be stored
	CreaContenedor(contenedor,run);			///< Creates the directory of the data and write a file with some parameters of the simulation in it. 
	
Float1D_MP meanDensity;			///< Instance of Float1D_MP to keep track of the mean density evolution
	InicializaFloat1D_MP(&meanDensity, T_max+10);	///<  initialization
	
Float1D_MP meanSize;			///< Instance of Float1D_MP to keep track of the mean size evolution
	InicializaFloat1D_MP(&meanSize, T_max+10);		///< initialization
	
float time_map[T_max+10];	///< To map monte-carlo sweep with fisical time
time_map[0]=0.0;		///< initialization

Rate_log Grate[6];			///< Structure that can be passed to the monte carlo sweep to keep track of instantaneous rates. (used to mesoure whatever you need)
InitRate_log(&Grate[5],10);	///< Initialization
InitRate_log(&Grate[0],50);	///< Initialization

FILE *file;		///< handle the files to write data to disk

/******************************************************/

	#pragma omp parallel			///< begin of parallel code. This code will be done for each thread 
	{
		int num_threads = omp_get_num_threads();	///< number of threads running
		int id = omp_get_thread_num();				///< id of 'this' thread
		int MaxPar = NoEnsambles/num_threads;		///< Number of simulation that each thread will handle
		#pragma omp master							///< The remainer simulations will be handle by the master thread
		{
				MaxPar+= NoEnsambles - MaxPar * num_threads;	 ///< Adding the remainer simulations to the master thread
		}
		estado e[MaxPar];						///< Struct "estado" stores the state of each simulation (positions, sizes, time)
		int Par;				///< Dummy variable used for looping.
		for(Par=0;Par<MaxPar;Par++)				///< for each simulated field
		{
			AlojaMemoria(&e[Par], NDX, NDY);	///< allocate memory for the field (state of the simualtion: positions, sizes, time)
			ResetEstado(&e[Par]);		///< initialization 
			e[Par].units=run.grid_units;	///< each simulation could have it's own units (not really used)
			e[Par].size_units=run.size_units; ///< each simulation could have it's own units (not really used)
		}
		init_JKISS(); ///< Initialization of the random number generator "JKISS"
			
		int i,j;		///< dummy variables used for looping
		
 /** INITIAL POSITIONS of individuals */
			for(Par=0;Par<MaxPar;Par++)			///< For each fiel (simulation) generate the INITIAL STATE
			{
		/** Random initial position with specific mean and minimun distance **/
		//	GeneraEstadoAleatorioTamano(&e[Par], 1.0/(float)run.initialMeanDistance , indv);	///< Random initial specifing the mean 
		//	FilterMinDistance(&e[Par], run.initialMinSeparation);	///< filter of minimun distance
		/**********************************************************/
		/** Grid initial positions **/
				//for(i=separation;i<NDX;i+=30*separation)
				//{
					//for(j=separation;j<NDY;j+=30*separation)
					//{
						//if(i+separation < NDX && j<NDY){
							//InsertIndividualAt(&e[Par],i,j,indv,1);
							//InsertIndividualAt(&e[Par],i+separation,j,indv,1);
						//}
					//}
				//}
		/*********************************************************/		
		/** Random initial positions specifing number of individuals */
				while(e[Par].ON < 55)
				{
					i=I_JKISS(1,NDX);
					j=I_JKISS(1,NDY);
					InsertIndividualAt(&e[Par],i,j,indv,0);
				}
		/**************************************************************/
			
		setMaxMetabolic(&e[Par],&modelo);		/// internal use. sets the frecuency of the poisson events in the simulation. 
		}

		Float1D_MP TamDist;				///< Temporal array to store size distribution to take advantage of each core cache. 
		InicializaFloat1D_MP(&TamDist, T_max+100);
			
	////////////////////////////// Monte CARLO :
	/** Initialization of temporal internal variables of each thread to improve memory efficiency **/
		char distT[50];
		Rate_log rate[6];
		int MaximoTamanoIni = 100;
		InitRate_log(&rate[0],MaximoTamanoIni);
		InitRate_log(&rate[1],MaximoTamanoIni);
		InitRate_log(&rate[2],MaximoTamanoIni);
		InitRate_log(&rate[3],MaximoTamanoIni);
		InitRate_log(&rate[4],MaximoTamanoIni);
		InitRate_log(&rate[5],rate[2].i_max+100);		
		indv.species=1;

	/** Loop for each time steep **/
		for(i=1;i<=T_max;i++)
		{			
			for(Par=0;Par<MaxPar;Par++) ///< loop for each simulation
			{
				BarrMCcRyCampTamano(&e[Par], run.Model.resource_rate, &modelo, rate); ///< A monte carlo sweep			
				ActualizaDistTamano_MP(&e[Par], &TamDist, 'A');				///< Store the size distribution in TamDist
				#pragma omp atomic
				meanDensity.array[e[Par].T]+=e[Par].ON;  		///< Store the mean density of all the simulations
			}			
			#pragma omp barrier
			#pragma omp atomic
			meanSize.array[TamDist.T]+=TamDist.array[0]; 	///< Store the the mean size of all the simulations
			#pragma omp master
			{
				time_map[e[0].T]=e[0].Meta_T;		///< Store the map between a computer time steep and fisical time.
			}
			#pragma omp barrier
				if((i-(i/write_interval)*write_interval)==1)    ///< write to disk just at certain time steeps
				{
					#pragma omp single
					{
						ResetFloat1D_MP(&TamDist_1);
					}
					SumaFloat1D_MP(&TamDist,&TamDist_1);		///< Add the size distribution of each thread to a global size distribution
					#pragma omp barrier
					#pragma omp single
					{
						sprintf(distT,"DT_%d",i);			///< Name of the file where the size distribution will be written
						GuardaFloat1D_MP(contenedor,distT,&TamDist_1); ///< Write to disk the global size distribution (all threads and simulations)
					}
					// estate
					#pragma omp master
					{
						GuardaEstadoEn(contenedor, &e[0]);		///< Write to disk the positions and sizes of the individuals of a single field
					}
					//
					init_JKISS();		///< Re-initialize the random number generation with a true random seed
				}

		}
	/**********************************************/

	#pragma omp single			///< A single thread should initialize the instance of the global rate 
	{
		FreeRate_log(&Grate[0]);
		InitRate_log(&Grate[0],rate[0].i_max);
		InitRate_log(&Grate[1],rate[1].i_max);
		InitRate_log(&Grate[2],rate[2].i_max);
		InitRate_log(&Grate[3],rate[3].i_max);
		InitRate_log(&Grate[4],rate[4].i_max);
	}

	/** Sum each thread rate to a global rate **/
	SumRate_log(&rate[0], &Grate[0]); 
	SumRate_log(&rate[1], &Grate[1]);
	SumRate_log(&rate[2], &Grate[2]);
	SumRate_log(&rate[3], &Grate[3]);
	SumRate_log(&rate[4], &Grate[4]);
	/****************************************/

	/** Free memory of the thread rate structure **/
	FreeRate_log(&rate[0]);
	FreeRate_log(&rate[1]);
	FreeRate_log(&rate[2]);
	FreeRate_log(&rate[3]);
	FreeRate_log(&rate[4]);
	FreeRate_log(&rate[5]);
	/********************************************/
	
		for(Par=0;Par<MaxPar;Par++)
		{
			LiberaMemoria(&e[Par]);  ///< Free memory of the array state for each simulation
		}
		
	
	}
	/****************END OF PARALLEL CODE *************************************/

/// Write calculated Rates to disk
FILE *taza;
int r;
char distT[50];
		sprintf(distT,"%s/GrowthR",contenedor);  ///< Name of the file where the evolution of the growth rate will be written.
		taza=fopen(distT, "w");						
		for(r=1;r<=Grate[0].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[0].Growth[r]/((double)Grate[0].NoEnsambles[r])); ///< Write the growth rate to file
		}
		fclose(taza);
		
		sprintf(distT,"%s/DeadR",contenedor); ///< Name of the file where the evolution of the dead rate will be written.
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[1].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, Grate[1].Growth[r]/((double)Grate[1].NoEnsambles[r])); ///< Write the dead rate to file
		}
		fclose(taza);
		
		sprintf(distT,"%s/ResourceR",contenedor);  ///< Name of the file where the evolution of the resource intake will be written 
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[2].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[2].Growth[r]/((float)Grate[2].NoEnsambles[r])); ///< Write the resource intake to file
		}
		fclose(taza);
		
		sprintf(distT,"%s/MetabolicR",contenedor);  
		taza=fopen(distT, "w");
		for(r=1;r<=Grate[3].i_max;r++){
			fprintf(taza,"%f %f\n",((float)r)*delta_s, delta_s*Grate[3].Growth[r]/((float)Grate[3].NoEnsambles[r]));
		}
		fclose(taza);
	
/** Free memory of Grate instance **/
FreeRate_log(&Grate[0]);	
FreeRate_log(&Grate[1]);
FreeRate_log(&Grate[2]);
FreeRate_log(&Grate[3]);
FreeRate_log(&Grate[5]);
/********************************/

LiberaMemoriaFloat1D_MP(&TamDist_1); ///< Free memory of TamDist_1 instance of Float1D

char thinning[50];
sprintf(thinning,"%s/thinning",contenedor); 							///< Open file to write some other rates
file=fopen(thinning, "w");												///
fputs("# T meanSize meanDensity fisicalTime totalResourceR\n",file);	///

for(r=1;r<=T_max;r++){			///< Write other rates to disk
fprintf(file,"%d %f %f %f %f\n",r, delta_s*(meanSize.array[r]/NoEnsambles), meanDensity.array[r]/NoEnsambles, time_map[r], Grate[4].Growth[r]/(coagulation_units*(float)Grate[4].NoEnsambles[r]) );
}
fclose(file);

FreeRate_log(&Grate[4]);			///< Free memory for Grate[4]
free(run.Model.meta_needs);			///< Free memory for meta_needs
return;
}
