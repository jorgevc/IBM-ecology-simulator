/*
Copyright 2012 Jorge Velazquez
*/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <omp.h> 
#include "libPP_5.0.h"
#include "EntSalArb_MP.h"
#include "ControlDinamico.h"

interfase *flags;
interfase flags_desconectado;
int shmid;
void *mem = NULL;

int tiempos[MAXNUCLEOS][2];

int* activa_escucha(int clave)
{
	
	if(clave==-1)
	{
		clave=getuid();
	}
		
	if((shmid = shmget((key_t)clave, sizeof(interfase), 0660 | IPC_CREAT)) == -1)
	{
		fprintf(stderr, "shmget fallo y no se pudo crear interfase para comandos\n");	
		flags_desconectado.guarda_estado = 0;
		flags_desconectado.guarda_RhoVsT = 0;
		flags_desconectado.servidor_ocupado = 0;
		flags_desconectado.mostrar_tiempos = 0;
		flags=&flags_desconectado;
	}else
	{
		if((mem = shmat(shmid, NULL, 0))==(void *)-1)
		{
			fprintf(stderr, "shmat fallo\n");
			flags_desconectado.guarda_estado = 0;
			flags_desconectado.guarda_RhoVsT = 0;
			flags_desconectado.servidor_ocupado = 0;
			flags_desconectado.mostrar_tiempos = 0;
			flags=&flags_desconectado;
		}else
		{
			flags = (interfase *)mem;
			flags->guarda_estado = 0;
			flags->guarda_RhoVsT = 0;
			flags->servidor_ocupado = 0;
			flags->mostrar_tiempos = 0;
			printf("%d\n",clave);
		}
	}
	int i;
	for(i=0;i<MAXNUCLEOS;i++)
	{
		tiempos[i][0]=0;
		tiempos[i][1]=0;
		flags->tiempos[i][0]=0;
		flags->tiempos[i][1]=-1;
	}
	
	return &(flags->servidor_ocupado);
}

void cierra_escucha(void)
{
	flags->servidor_ocupado=0;
	
	 if (shmdt(mem) == -1) {
        fprintf(stderr, "shmdt detach fallo\n");
    }

    if (shmctl(shmid, IPC_RMID, 0) == -1) {
        fprintf(stderr, "shmctl(IPC_RMID) fallo: eliminar memoria compartida\n");
    }

return;	
}

void Servicio(int T,char *contenedor,Float2D_MP *MP_RhoVsT, Float2D_MP *MP_RhoVsT_1, Dist_MP *MP_RhoDist, Dist_MP *MP_RhoDist_1, estado *e, int MaxPar)
{
	volatile static int maximo_listo=0;
	volatile static int tiempo_maximo=-2;
	int id = omp_get_thread_num();
	int num_threads = omp_get_num_threads();
	int tiempos_listos=0;
	int Par,n;

	if(flags->mostrar_tiempos==1)
	{
		printf("Thread id: %d en T=%d\n",id,T);
		flags->tiempos[id][0]=T;
		flags->tiempos[id][1]=1;
		if(id==0)
		{
			while(tiempos_listos<num_threads)
			{
				tiempos_listos=0;
				for(n=0;n<num_threads;n++)
				{
					tiempos_listos+=flags->tiempos[n][1];
				}
				printf("esperando a que los otros threats por su tiempo ...\n");
				sleep(1);				
			}
			flags->mostrar_tiempos=0;
			flags->servidor_ocupado=0;
		}
		
		while(flags->mostrar_tiempos==1)
		{
			sleep(1);
		}
		
		flags->tiempos[id][1]=0;
	}
	
	if(flags->guarda_RhoVsT == 1 || flags->guarda_estado == 1)
	{
		if(tiempos[id][1]==0)
		{
			tiempos[id][0]=T;
			tiempos[id][1]=1;
		}
		
		#pragma omp flush(maximo_listo,tiempo_maximo)
		
		if(id==0 && maximo_listo!=1)
		{
			while(tiempos_listos<num_threads)
			{
				tiempos_listos = 0;
				for(n=0;n<num_threads;n++)
				{
					tiempos_listos+=tiempos[n][1];
				}
				printf("esperando a que los otros threats a darme su tiempo ...\n");
				sleep(1);
			}
			
			for(n=0;n<num_threads;n++)
			{
				if(tiempo_maximo<tiempos[n][0])
				{
					tiempo_maximo=tiempos[n][0];
				}
			}
		maximo_listo=1;	
		#pragma omp flush(maximo_listo,tiempo_maximo)
		printf("Tiempo maximo listo en T=%d\n",tiempo_maximo);
		}
		
		while(maximo_listo!=1)
		{
			printf("therad:%d, esperando a que este listo el tiempo de la barrera \n",id);
			sleep(1);	
		}
		
			if(T==(tiempo_maximo + 1))
			{
				#pragma omp barrier
				if(flags->guarda_RhoVsT==1)
				{
					SumaFloat2D_MP(MP_RhoVsT, MP_RhoVsT_1);
					for(Par=0;Par<MaxPar;Par++)
					{
						ActualizaRhoVsT_MP(&e[Par],NULL,MP_RhoDist);
					}
					SumaDist_MP(MP_RhoDist, MP_RhoDist_1);
					#pragma omp barrier
					#pragma omp single
					{
						puts("Guardando informacion...\n");
						GuardaRhoVsT_MP(contenedor,MP_RhoVsT_1,MP_RhoDist_1);
						GuardaTiposEn_MP(contenedor,MP_RhoVsT_1,T);
						flags->guarda_RhoVsT=0;
						sprintf(flags->mensaje,"RhoVsT Guardado a T=%d\n",T);
						ResetFloat2D_MP(MP_RhoVsT_1);
						ResetDist_MP(MP_RhoDist_1);
					}
				}
				if(flags->guarda_estado==1)
				{
					for(Par=0;Par<MaxPar;Par++)
					{
						GuardaEstadoEn_MP(contenedor,&e[Par],id,Par);
					}
					#pragma omp barrier
					if(id==0)
					{
						flags->guarda_estado=0;
						sprintf(flags->mensaje,"Estado Guardado a T=%d\n",T);
					}
				}
				
				tiempos[id][1]=0;
				#pragma omp barrier
				if(id==0)
				{
					maximo_listo=0;
					flags->servidor_ocupado=0;
				}
				#pragma omp barrier
			}
			
	}
		
	
	return;
}

void SalidaCD(int *i,int T_max)
{
	if(flags->servidor_ocupado==1)
	{
		if(flags->salida_loop==1)
		{
			*i=T_max+1;
		}
		if(flags->salida_violenta==1)
		{
			#pragma omp single
			{
				cierra_escucha();
				exit(0);
			}
		}
	}
	return;
}


#include <signal.h>

/* This flag controls termination of the main loop.  */
volatile sig_atomic_t keep_going = 1;
volatile sig_atomic_t guarda_est = 0;
#pragma omp threadprivate(guarda_est)

/* The signal handler just clears the flag and re-enables itself.  */
void catch_alarm (int sig)
{
  keep_going = 0;
  signal (sig, catch_alarm);
}

void guarda_estado(int sig)
{
	guarda_est=1;
	signal(sig, guarda_estado);
}

/*
#pragma omp master
	{
		signal (SIGALRM, catch_alarm);
	   Set an alarm to go off in a little while.  
		alarm (TIEMPO_ESPERA);
	}

#pragma omp barrier
			if(keep_going==0)
			{
				#pragma omp master
				{
					keep_going=1;
					alarm (TIEMPO_ESPERA);
				}
				for(Par=0;Par<MaxPar;Par++)
				{
					GuardaEstadoEn_MP(contenedor,&e[Par],id,Par);
				}
			}
*/
