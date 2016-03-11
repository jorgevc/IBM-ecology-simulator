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

#include "MLE.c"

main(){	
#ifdef VIRTUAL_GRID
printf("VIRUAL_GRID=TRUE\n");
#endif
#ifdef SOI
printf("SOI=TRUE\n");
#endif

///////////////////////////Inicializa parametros de la simulacion
runDescriptor run;
run.X=200;  	//unidades "fisicas"
run.Y=run.X;
run.grid_units=1.0; //factor de conversion de unidades "fisicas" a lado de celda
run.size_units=1.0; //numero de unidades en computo que hacen una unidad "fisica" de tamano (divisor de conversion)

//run.T_max=(run.grid_units*run.grid_units)*802;
//run.T_max=(run.grid_units*run.grid_units)*50000;
//run.T_max=(run.grid_units*run.grid_units)*20000;
run.T_max=3000;
run.NoEnsambles=4;
int const write_interval=2000;

//in this case the Area_units and Length_units are units of the imaginary resources grid.
run.Model.ResourcesScale = 1.0; //conversion factor to the imaginary resources grid or SOI resolution
//run.Model.competitionAsymetry = (8.0/3.0)*1.6212;
run.Model.competitionAsymetry = 80.0;
double const Area_units=(run.grid_units*run.Model.ResourcesScale)*(run.grid_units*run.Model.ResourcesScale);
double const Length_units=run.grid_units*run.Model.ResourcesScale;

double coagulation_units = 1.0;
#ifdef SOI
coagulation_units = 52.0;
run.Model.coagulation_factor = coagulation_units*1.0;
#else
run.Model.coagulation_factor = coagulation_units*1.0*3.1416;
#endif
//run.Model.coagulation_factor=((Area_units/pow(run.size_units,run.Model.coagulation_exp))*1.0);

run.Model.coagulation_radio_exp=1.0; //cambiar tambien en model.c  :0.375
run.Model.coagulation_exp=2*run.Model.coagulation_radio_exp;  // cambiar tambien en model.c
run.Model.Cr=1.0;
run.Model.coagulation_radio_factor=1.0;
run.Model.metabolic_exp=2.6666; //cambiar tambien en model.c
float dgrmax=48.6;
run.Model.Cm=(3.1416*pow(2.0,2.0/3.0)*pow(run.Model.Cr,8.0/3.0))/(3.0*pow(dgrmax,2.0/3.0));
float dmax=2.0*pow(3.1416,3.0/2.0)*pow(run.Model.Cr,4.0)/pow(run.Model.Cm,3.0/2.0);
run.Model.metabolic_factor=(run.Model.Cm)/pow(2.0*run.Model.Cr,8.0/3.0);
run.Model.health_factor=0.0; //usandose lineal proporcional al tamano (adimensional) fraccion de biomasa que puede "danarse" antes de enfermar. 
// con 2.0 salen buenos resultados para U shape de dead rate.
run.Model.Cg=4.0;
run.Model.growth_constant=(3.0*pow(run.Model.Cr,8.0/3.0))/(pow(2.0,1.0/3.0)*run.Model.Cg); //  needed resources per unit size increse. 
//run.initialMeanDistance=run.grid_units*12;
//run.initialMinSeparation=1000;
run.scaleFactor=10.00;
//run.Model.meta_needs=SetMetaNeeds(run.Model, Area_units);
//run.Model.R=SetR(run.Model, Length_units);
//run.Model.M=SetM(run.Model, Area_units);

#ifdef HEALTH_TRACK	
run.Model.min_health=0;
#else
run.Model.min_health=0;
#endif


//run.Model.birth_rate=0.5;
run.Model.birth_rate=0.0;
run.Model.RadioBirth=15;
run.Model.dead_rate=0.0;
run.Model.intra_coagulation=0.0;

float like,a,b,c,d,e,f,aceptedLike,r,h,dgrmaxA, tmpCr;
int i,g,n;
i=0;
char contenedor[150];
	sprintf(contenedor,"DATOS_MCMC/Ene28_MinVariance_climbing93");
	CreaContenedor(contenedor,run);

	
char contenedorLec[150];	
sprintf(contenedorLec,"OBSERVACIONES/COOMESDATA/commes.csv");
DataSet* datos[250];
for(n=0;n<250;n++)
{
	datos[n] = loadSizes(contenedorLec, n+1);
//	printf("\nNo de datos leidos %d\n", datos[n]->noPoints);	
}
run.initialPopulation = datos[0]->noPoints;
		aceptedLike=0.0;
		dgrmaxA=dgrmax;
		a=run.Model.Cm;
		b=run.Model.coagulation_radio_factor;
		c=run.Model.competitionAsymetry;
		d=run.Model.resource_rate;
		e=run.Model.Cg;
		f=run.Model.health_factor;
		g=run.initialPopulation;
		h=run.scaleFactor;
		tmpCr=run.Model.Cr;

char chain[150];
sprintf(chain,"%s/chain",contenedor);
FILE *file;
file=fopen(chain, "r");
if(file == NULL)
{
file=fopen(chain, "a");	
//fputs("# No,likelihood,Cm,Coagulation_radio_factor,Competition_asymetry,Resource_rate,Cg,Health_factor,Initial_population,Scale_factor, time84, time93\n",file);
fputs("# No, likelihood, Cm, CompetitionAsymetry, Scale, Cr ,Cg, time\n",file);
fclose(file);
i=0;
aceptedLike=0.0;
}else{
	char *buffer;
	size_t tam_buffer = 300*sizeof(char);
	buffer = (char *) malloc (tam_buffer);	
	if(0==fseek(file, -600, SEEK_END))
	{
		getline(&buffer, &tam_buffer, file);
		int args_assigned = 0;
		while(getline(&buffer, &tam_buffer, file)!=-1)
		{
			if(strchr(buffer, '#')==NULL)
			{
				args_assigned = sscanf(buffer, "%d %f %f %f %f %f %f %f %d %f %*d %*d",&i,&aceptedLike,&a,&b,&c,&d,&e,&f,&g,&h);
				if(args_assigned == 10)
				{
					run.Model.Cm=a;
					run.Model.coagulation_radio_factor=b;
					run.Model.competitionAsymetry=c;
					run.Model.resource_rate=d;
					run.Model.Cg=e;
					run.Model.health_factor=f;
					run.initialPopulation=g;
					run.scaleFactor=h;
				}
			}
		}
	}
	free(buffer);
	fclose(file);
}

		
		
	r=1.0;
int time84,time93;
Pair partial,cum;
while(1)
{
	i++;
		do
		{
			run.Model.Cr=tmpCr + (F_JKISS() - 0.5);
		}while(run.Model.Cr<=0.0);
		do
		{
			//dgrmax=dgrmaxA - 20.0*(F_JKISS() - 0.5);
			run.Model.Cm= a + (F_JKISS() - 0.5);
			run.Model.metabolic_factor=(run.Model.Cm)/pow(2.0*run.Model.Cr,8.0/3.0);	
		}while(run.Model.metabolic_factor<=0.0);
	//	do
	//	{
	//		run.Model.coagulation_radio_factor = b + 0.01*2.0*(F_JKISS() - 0.5);  //checar cuando cambie size units o coagulation_exp
	//	}while(run.Model.coagulation_radio_factor<=0.0);
		do
		{
			run.Model.competitionAsymetry = c - 20.0*(F_JKISS() - 0.5);
		}while(run.Model.competitionAsymetry < 0.01);
		//do
		//{
			//run.initialPopulation=g + 10*2.0*(F_JKISS() - 0.5);
		//}while((run.initialPopulation + (run.initialPopulation/10)) < datos->noPoints );
		//do
	//	{
	//		run.scaleFactor= h + 0.1*2.0*(F_JKISS() - 0.5);
	//	}while(run.scaleFactor < 0.01);
		do
		{
			run.Model.Cg = e + (2.0*(F_JKISS() - 0.5));
			run.Model.growth_constant=(3.0*pow(run.Model.Cr,8.0/3.0))/(pow(2.0,1.0/3.0)*run.Model.Cg);
		}while(run.Model.growth_constant <= 0.0);
	//	do
	//	{
	//		run.Model.health_factor = f + 0.2*2.0*(F_JKISS() - 0.5);
	//	}while(run.Model.health_factor < 0.01);
		//do
		//{
			//run.Model.resource_rate = d + 0.3*2.0*(F_JKISS() - 0.5);
		//}while(run.Model.resource_rate < 0.1);
		
	printf("cr=%f, cm=%f, cg=%f \n",run.Model.Cr,run.Model.Cm,run.Model.Cg);	
		
	cum.A=0;
	cum.B=0;	
	for(n=0;n<250;n++)
	{	
		run.initialPopulation = datos[n]->noPoints;
		partial = MLE(datos[n], run, &time84, &time93);
		cum.A += partial.A;
		cum.B += partial.B; 
	}
	
	like = ((float)cum.A/(float)cum.B);
	
	if(aceptedLike < 0){
		r=aceptedLike/like;
	}
	if(r>=1.0)
	{
		r = 1.0;
	}else{
		r*=0.2*r;
	}
	if(F_JKISS() < r)
	{
		aceptedLike=like;
		a=run.Model.Cm;
		b=run.Model.coagulation_radio_factor;
		c=run.Model.competitionAsymetry;
		d=run.Model.resource_rate;
		e=run.Model.Cg;
		f=run.Model.health_factor;
		g=run.initialPopulation;
		h=run.scaleFactor;
		dgrmaxA=dgrmax;
		tmpCr=run.Model.Cr;
		file=fopen(chain,"a");
		//fprintf(file,"%d %f %f %f %f %f %f %f %d %f %d %d\n",i,aceptedLike,a,b,c,d,e,f,g,h,time84, time93);
		fprintf(file,"%d %f %f %f %f %f %f %d\n",i,aceptedLike,a,c,h,tmpCr,e,time84);
		fclose(file);
		printf("iteracion %d Aceptada! like:%f\n",i,aceptedLike);
	}else{
		printf("iteracion %d Rechazada! like:%f\n",i,like);
	}			
}

for(n=0;n<250;n++)
{
	FreeDataSet(datos[n]);
}
return;
}
