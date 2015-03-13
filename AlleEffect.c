/*
Copyright 2015 Jorge Velazquez
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
#include <math.h>
#include <time.h>

//#include "libPP_5.1.h"
#include "MC_alle_effect.h"
#include "EntSalArb_MP.h"
#include "GNA.h"
#include "conn_mysql.h"


main(){	
	
///////////////////////////Inicializa parametros de la simulacion
int NDX=150;
int NDY=NDX;
int T_max = 100;
int NoEnsambles=20;

int CantidadEspecies=2;

alle_env *env = (alle_env *)malloc(1 * sizeof(alle_env));
env->param = (especie *)malloc((CantidadEspecies + 1) * sizeof(especie));
env->param[1].Birth= 1.0;
env->param[1].Coagulation = 1.0; //Brown usa: 0.00002; 
env->param[1].CoagulationIntra= 0.2; //Modelo J-C 0.0008
env->param[1].Dead= 0.4;
env->param[1].RadioBirth= 1;
env->param[1].RadioCoa= 10;
env->param[1].RadioCoaIntra= 10;  //Modelo Heteromyopia 20


env->param[2].Birth= 1.0;
env->param[2].Coagulation= 0.2; //Brown usa: 0.00002;
env->param[2].CoagulationIntra=0.2;
env->param[2].Dead = 0.4;
env->param[2].RadioBirth= 1;
env->param[2].RadioCoa= 10;
env->param[2].RadioCoaIntra= 10;


omp_set_num_threads(4);

//date time of simulation
time_t raw_now = time(NULL);
struct tm * now;
now = localtime ( &raw_now );
char sim_time[50];
sprintf(sim_time,"'%d-%d-%d %d:%d:%d'",now->tm_year + 1900, now->tm_mon + 1, now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);


////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:

Float2D_MP MP_RhoVsT_1;		
		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, CantidadEspecies, 0);
		
Float1D_MP MP_Correlacion_1G;
		InicializaFloat1D_MP(&MP_Correlacion_1G, NDX);
		
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS	
			
			///////////////////////////////////// INICIA PARALLEL

			#pragma omp parallel			///////Estado INICIAL:
			{
				init_JKISS(); //Inicializa la semilla de cada proceso.
				
				int num_threads = omp_get_num_threads();
				int id = omp_get_thread_num();
				int MaxPar = NoEnsambles/num_threads;
				#pragma omp master
				{
						MaxPar+= NoEnsambles - MaxPar * num_threads;	 
				}
				estado e[MaxPar];
				
				int Par;
				for(Par=0;Par<MaxPar;Par++)
				{
					AlojaMemoria(&e[Par], NDX, NDY);  
					ResetEstado(&e[Par]); 	
				}

				
					for(Par=0;Par<MaxPar;Par++)
					{
					//InsertaIndividuosAleatorio(&e[Par],100,NoEspecie);
					GeneraEstadoAleatorio(&e[Par], 0.2, 1);
					GeneraEstadoAleatorio(&e[Par], 0.2, 2);
					}


				
			///////////////////////////////////////Termina Estado INICIAL
			////////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

			Float2D_MP MP_RhoVsT;	
				InicializaFloat2D_MP(&MP_RhoVsT, T_max, CantidadEspecies, MaxPar);
						
			Float2D_MP MP_Corr2D_1;
			InicializaFloat2D_MP(&MP_Corr2D_1, NDX, NDY, 0);
			
			Float1D_MP MP_Correlacion_1;
			InicializaFloat1D_MP(&MP_Correlacion_1, NDX);
								
			///////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO
					
			////////////////////////////////Barrido Monte CARLO:
				int i;
				for(i=0;i<T_max;i++)
				{
					for(Par=0;Par<MaxPar;Par++)
					{
						MC_sweep_alle(&e[Par], env);
					//necesario?	ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);	
					}	
				
						if((i-(i/500)*500)==499)    //Inicializa cada 500 pasos
						{
							init_JKISS();
						}
		
						
						///////////////////////Evolucion de Correlacion
						
						//if((i-(i/30)*30)==0)    // cada 30 pasos
						//{
							//ResetFloat1D_MP(&MP_Correlacion_12G);
							//ResetFloat1D_MP(&MP_Correlacion_12);
							//ResetFloat2D_MP(&MP_Corr2D_12);	
							//CFFT_Tipos_MP(e, MaxPar, &MP_Corr2D_12, 1, 2);
							
							//if(MP_Corr2D_12.NoEnsambles > 0)
							//{
								//CompactaCorrelacion(&MP_Corr2D_12, &MP_Correlacion_12);
								//SumaFloat1D_MP(&MP_Correlacion_12,&MP_Correlacion_12G);
							//}
							//#pragma omp barrier
							//#pragma omp master
							//{
								//if(i>400)
								//{
									//GuardaCorrelacion_MP(contenedorCompleto, "1-2" , &MP_Correlacion_12G);
								//}
							//}
						//}
				}
				
				
			////////////////////////////////Termina Monte CARLO
				for(Par=0;Par<MaxPar;Par++)
				{
					ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
				}
				
			#pragma omp barrier
			#pragma omp single
			{
				ResetFloat2D_MP(&MP_RhoVsT_1);
			}
			
				SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
				
			//	#pragma omp master
				//{
				//	PD_GuardaEstadoEn_MP(contenedorCompleto, e, id, 1);
				//}
			////////Correlacion
			
			
			//ResetFloat2D_MP(&MP_Corr2D_1);	
			//ResetFloat1D_MP(&MP_Correlacion_1);
									
			//ResetFloat1D_MP(&MP_Correlacion_1G);
			
						//CFFT_Tipos_MP(e, MaxPar, &MP_Corr2D_1, 1, 1);
						
						//if(MP_Corr2D_1.NoEnsambles > 2)
						//{
							//CompactaCorrelacion(&MP_Corr2D_1, &MP_Correlacion_1);
							//SumaFloat1D_MP(&MP_Correlacion_1,&MP_Correlacion_1G);
						//}
						
						//Libera Memoria
						for(Par=0;Par<MaxPar;Par++)
						{
							LiberaMemoria(&e[Par]);
						}
						LiberaMemoriaFloat2D_MP(&MP_RhoVsT);	
						
						//LiberaMemoriaFloat2D_MP(&MP_Corr2D_1);
						
						//LiberaMemoriaFloat1D_MP(&MP_Correlacion_1);

			}	////////////////////////////////////////////////////////////////////TERMINA PARALLEL
	
			//// Guarda parametros en MySql	y crea CONTENEDOR
	char contenedor[300];
	sprintf(contenedor,"test");
	CreaContenedor(contenedor);	
		
	//char values[300];
	//sprintf(values,"NULL,%s,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%d,%d,%d,%d,'%s',0",sim_time,
	//...
	//NDX,
	//NDY,
	//T_max,
	//NoEnsambles,
	//contenedor
	//);
	
	//printf("MySQL client version: %s\n", mysql_get_client_info());
	//MYSQL *con = connect_db("localhost","usr", "password", "db");
	//int inserted_id = insert_into_db(con, "table",values);
	//mysql_close(con);
	char contenedorCompleto[200];
	//sprintf(contenedorCompleto,"%s/%d",contenedor,inserted_id);
	sprintf(contenedorCompleto,"%s/%d",contenedor,1);
	CreaContenedor(contenedorCompleto);
	
	
	//store_density_evolution(contenedorCompleto,&MP_RhoVsT_1, 0, FisicalTime,R_dyn,Temp_dyn);
	//GuardaCorrelacion_MP(contenedorCompleto, "1-1" , &MP_Correlacion_1G);			
	GuardaRhoVsT_MP(contenedorCompleto,&MP_RhoVsT_1,NULL);	
	
	//FILE *aA;
	//char archivo[200];
	//sprintf(archivo,"Graficas/simulation.tex");	
	//aA=fopen(archivo, "w");
	//fprintf(aA,"\\newcommand{\\data}{../%s/density_evolution}\n\\newcommand{\\plotTitle}{id=%d  }",contenedorCompleto,inserted_id);
	//fclose(aA);
	//sprintf(archivo,"Graficas/include_make");
	//aA=fopen(archivo, "w");
	//fprintf(aA,"id = %d\n",inserted_id);
	//fclose(aA);
	//system("cd Graficas; make figure");
	//////
	
 
	LiberaMemoriaFloat2D_MP(&MP_RhoVsT_1);

						
return;
}
