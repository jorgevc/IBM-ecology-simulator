/*
Copyright 2012 Jorge Velazquez
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "libPP_5.0.h"
#include "EntSalArb_MP.h"
#include "GNA.h"

main(){	
	
///////////////////////////Inicializa parametros de la simulacion
int NDX=150;
int NDY=NDX;
int T_max = 10000;
int NoEnsambles=20;

int CantidadEspecies=2;


float Birth1= 1.0;
float Coagulation1; //Brown usa: 0.00002; 
float CoaIntra1= 0.2; //Modelo J-C 0.0008
float Dead1= 0.4;
int RadioBirth1= 1;
int RadioCoa1= 10;
int RadioCoaIntra1= 10;  //Modelo Heteromyopia 20


float Birth2= 1.0;
float Coagulation2= 0.2; //Brown usa: 0.00002;
float CoaIntra2=0.2;
float Dead2;
int RadioBirth2= 1;
int RadioCoa2= 10;
int RadioCoaIntra2= 10;


omp_set_num_threads(4);
////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:
int T_st;
float Rho_A;

Float2D_MP MP_RhoVsT_1;
int NoEspecies=CantidadEspecies;		
		InicializaFloat2D_MP(&MP_RhoVsT_1, T_max, CantidadEspecies, 0);
		
Float1D_MP MP_Correlacion_1G;
		InicializaFloat1D_MP(&MP_Correlacion_1G, NDX);
		
Float1D_MP MP_Correlacion_2G;
		InicializaFloat1D_MP(&MP_Correlacion_2G, NDX);

Float1D_MP MP_Correlacion_12G;
		InicializaFloat1D_MP(&MP_Correlacion_12G, NDX);
		
char contenedor[300];
	
FILE *pD;
char NombrePD[200];
char contenedorCompleto[320];
/////////////////////////////////////Termina Prepara CONTENEDOR para escribir DATOS

RadioCoa1=10;
RadioCoaIntra1=10;

RadioCoa2=10;
RadioCoaIntra2=10;
	
	sprintf(contenedor,"DATOS2/PD(0.2:0.2)_(B=1:1,D=%1.2f:d2,C=c21:%1.2f,CI=%1.2f:%1.2f,RB=1:1,RC=%d:%d,RCI=%d:%d)_(NDX=%d,Tmax=%d)",Dead1,Coagulation2,CoaIntra1,CoaIntra2,RadioCoa1,RadioCoa2,RadioCoaIntra1,RadioCoaIntra2,NDX,T_max);
	CreaContenedor(contenedor);
	
	sprintf(NombrePD,"%s/PD",contenedor);
	pD=fopen(NombrePD, "a");
	fprintf(pD,"#d2 c21 Rho_1 Rho_2 time\n"); 
	fclose(pD);
	
  for(Dead2=0.25;Dead2<0.52;Dead2+=0.01)
  {
	  		//pD=fopen(NombrePD, "a");
			//fprintf(pD,"\n"); 
			//fclose(pD);
	for(Coagulation1=0.0;Coagulation1<0.52;Coagulation1+=0.01)
	{			
			
			SetSpecie2(1, Birth1, Coagulation1, CoaIntra1, Dead1, RadioBirth1, RadioCoa1, RadioCoaIntra1);
			SetSpecie2(2, Birth2, Coagulation2, CoaIntra2, Dead2, RadioBirth2, RadioCoa2, RadioCoaIntra2);
			
			//sprintf(contenedorCompleto,"%s/(d2=%1.2f,c21=%1.2f)*",contenedor,Dead2,Coagulation1);
			//CreaContenedor(contenedorCompleto);
	
			ResetFloat2D_MP(&MP_RhoVsT_1);
			
			
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
				
				//MaxPar=CargaEstado_MP(contenedorLec,"T_001",e,NDX,NDY,id,MaxPar);
				
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
						InicializaFloat2D_MP(&MP_RhoVsT, T_max, NoEspecies, MaxPar);
						
			//Float2D_MP MP_Corr2D_1;
			//InicializaFloat2D_MP(&MP_Corr2D_1, NDX, NDY, 0);
			
			//Float2D_MP MP_Corr2D_2;
			//InicializaFloat2D_MP(&MP_Corr2D_2, NDX, NDY, 0);
			
			//Float2D_MP MP_Corr2D_12;
			//InicializaFloat2D_MP(&MP_Corr2D_12, NDX, NDY, 0);
			
			//Float1D_MP MP_Correlacion_1;
			//InicializaFloat1D_MP(&MP_Correlacion_1, NDX);
			
			//Float1D_MP MP_Correlacion_2;
			//InicializaFloat1D_MP(&MP_Correlacion_2, NDX);
			
			//Float1D_MP MP_Correlacion_12;
			//InicializaFloat1D_MP(&MP_Correlacion_12, NDX);
								
			///////////////////////////////////Termina prepara Contenedor MEMORIA de cada PROCESO

					////for(Par=0;Par<MaxPar;Par++)
					////{
						////ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);	
					////}
					
			////////////////////////////////Barrido Monte CARLO:
				int i;
				for(i=0;i<T_max;i++)
				{
					for(Par=0;Par<MaxPar;Par++)
					{
						BarrMCcRyCamp(&e[Par]);
					}
					
					if((i-(i/200)*200)==0)
					{
						for(Par=0;Par<MaxPar;Par++)
						{
							ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);	
						}
					}
					if((i-(i/200)*200)==20)
					{
						for(Par=0;Par<MaxPar;Par++)
						{
							ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
						}
						SumaFloat2D_MP(&MP_RhoVsT, &MP_RhoVsT_1);
						
						#pragma omp barrier
						if(fabs(MP_RhoVsT_1.array[i+1][1]-MP_RhoVsT_1.array[i-19][1])<0.001 || fabs(MP_RhoVsT_1.array[i+1][2]-MP_RhoVsT_1.array[i-19][2])<0.001)
						{
							#pragma omp master
							{
								T_st=i+1;
							}
							i=T_max+1;
						}
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
			if(i==T_max)
			{
				for(Par=0;Par<MaxPar;Par++)
				{
					ActualizaRhoVsT_MP(&e[Par],&MP_RhoVsT,NULL);
				}
				#pragma omp master
				{
					T_st=T_max;
				}
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
			//ResetFloat2D_MP(&MP_Corr2D_2);
			//ResetFloat2D_MP(&MP_Corr2D_12);	
			
			//ResetFloat1D_MP(&MP_Correlacion_1);
			//ResetFloat1D_MP(&MP_Correlacion_2);
			//ResetFloat1D_MP(&MP_Correlacion_12);
									
			//ResetFloat1D_MP(&MP_Correlacion_1G);
			//ResetFloat1D_MP(&MP_Correlacion_2G);
			//ResetFloat1D_MP(&MP_Correlacion_12G);
			
						//CFFT_Tipos_MP(e, MaxPar, &MP_Corr2D_1, 1, 1);
						//CFFT_Tipos_MP(e, MaxPar, &MP_Corr2D_2, 2, 2);
						//CFFT_Tipos_MP(e, MaxPar, &MP_Corr2D_12, 1, 2);
						
						//if(MP_Corr2D_1.NoEnsambles > 2)
						//{
							//CompactaCorrelacion(&MP_Corr2D_1, &MP_Correlacion_1);
							//SumaFloat1D_MP(&MP_Correlacion_1,&MP_Correlacion_1G);
						//}
						//if(MP_Corr2D_2.NoEnsambles > 2)
						//{
							//CompactaCorrelacion(&MP_Corr2D_2, &MP_Correlacion_2);
							//SumaFloat1D_MP(&MP_Correlacion_2,&MP_Correlacion_2G);
						//}
						//if(MP_Corr2D_12.NoEnsambles > 2)
						//{
							//CompactaCorrelacion(&MP_Corr2D_12, &MP_Correlacion_12);
							//SumaFloat1D_MP(&MP_Correlacion_12,&MP_Correlacion_12G);
						//}
						
						//Libera Memoria
						for(Par=0;Par<MaxPar;Par++)
						{
							LiberaMemoria(&e[Par]);
						}
						LiberaMemoriaFloat2D_MP(&MP_RhoVsT);	
						
						//LiberaMemoriaFloat2D_MP(&MP_Corr2D_1);
						//LiberaMemoriaFloat2D_MP(&MP_Corr2D_2);	
						//LiberaMemoriaFloat2D_MP(&MP_Corr2D_12);
						
						//LiberaMemoriaFloat1D_MP(&MP_Correlacion_1);
						//LiberaMemoriaFloat1D_MP(&MP_Correlacion_2);
						//LiberaMemoriaFloat1D_MP(&MP_Correlacion_12);

			}	////////////////////////////////////////////////////////////////////TERMINA PARALLEL

					//GuardaCorrelacion_MP(contenedorCompleto, "1-1" , &MP_Correlacion_1G);
					//GuardaCorrelacion_MP(contenedorCompleto, "2-2" , &MP_Correlacion_2G);
					//GuardaCorrelacion_MP(contenedorCompleto, "1-2" , &MP_Correlacion_12G);
			
			//GuardaRhoVsT_MP(contenedorCompleto,&MP_RhoVsT_1,NULL);	
				if(MP_RhoVsT_1.array[T_st][1]>0.0 || MP_RhoVsT_1.array[T_st][2]>0.0)
				{
					pD=fopen(NombrePD, "a");
					fprintf(pD,"%f %f %f %f %d\n",Dead2,Coagulation1, MP_RhoVsT_1.array[T_st][1]/(float)MP_RhoVsT_1.NoEnsambles, MP_RhoVsT_1.array[T_st][2]/(float)MP_RhoVsT_1.NoEnsambles,T_st); 
					fclose(pD);
				}
				
				if(fabs(MP_RhoVsT_1.array[T_st][1] - Rho_A) < 2.0*0.01*20.0)
				{
					Coagulation1+=0.01;
				} 
				Rho_A=MP_RhoVsT_1.array[T_st][1];
	
	}
  
  }
 
						LiberaMemoriaFloat2D_MP(&MP_RhoVsT_1);

						
return;
}
