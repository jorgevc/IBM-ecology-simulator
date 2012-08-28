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

#include "InterfaseControlDinamico.h"

#include <signal.h>

void GuardaEstado(void);


/* This flag controls termination of the main loop.  */
volatile sig_atomic_t keep_going = 1;
volatile sig_atomic_t FLguarda_estado = 0;

/* The signal handler just clears the flag and re-enables itself.  */
void catch_alarm (int sig)
{
	FLguarda_estado = 1;
		signal (sig, catch_alarm);
}

int main()
{
	int FLcrono_establecido = 0;
    int running = 1;
    int Tiempo=600;
    void *mem = NULL;
    interfase *flags;
    char buffer[100];
    int shmid,id;
    int clave=getuid();
    signal(SIGALRM, catch_alarm);  
    
	printf("Clave de conexion:\n");
    fgets(buffer, 100, stdin);
    sscanf(buffer, "%d", &clave);
    printf("%d\n",clave);

    if((shmid = shmget((key_t)clave, sizeof(interfase), 0666 | IPC_CREAT)) == -1)
    {
		fprintf(stderr, "shmget fallo\n");
        exit(EXIT_FAILURE);
	}

    if((mem = shmat(shmid, NULL, 0))==(void *)-1)
    {
        fprintf(stderr, "shmat fallo\n");
        exit(EXIT_FAILURE);
    }

	printf("Enlazado al proceso con clave:%d\n",clave);
	
    flags = (interfase *)mem;
    while(running==1) {
		 while(FLcrono_establecido == 1)
        {
			if(FLguarda_estado == 1)
			{
				FLguarda_estado = 0;
				while(flags->servidor_ocupado == 1) 
				{
					printf("servidor ocupado...\n");
					sleep(1);            
				}
				flags->servidor_ocupado = 1;
				flags->guarda_estado = 1;
				printf(flags->mensaje);
				alarm(Tiempo);
			}
			sleep(10);
		}
		
        while(flags->servidor_ocupado == 1) {
			printf("servidor ocupado...\n");
            sleep(1);            
        }
       
        printf("Menu: \n");
        printf("Guarda estado: guarda todos los ensambles del experimento\n");
        printf("Guarda RhoVsT: guarda RhoVsT y DisRho. NO GUARDA LOS ENSAMBLES!\n");
        printf("Mostrar tiempos: muestra el tiempo en que va cada proceso\n");
        printf("Salir del loop principal: lo que su nombre indica\n");
        printf("Matar al proceso: lo que su nombre indica\n");
        printf("Reconectarse: conectarse con otro proceso\n");
        printf("Crono:[tiempo en seg]: Establecer Crono\n");
        printf("Elimina Crono:lo que su nombre indica\n");
        printf("salir: Sale de esta interfase\n");
        printf("Escribe la orden:\n");
        
        fgets(buffer, 100, stdin);
        
         if (strncmp(buffer, "Guarda estado", 13) == 0) {
                flags->guarda_estado = 1;
                flags->servidor_ocupado = 1;
                printf(flags->mensaje);
        }
         if (strncmp(buffer, "Guarda RhoVsT", 13) == 0) {
                flags->guarda_RhoVsT = 1;
                flags->servidor_ocupado = 1;
                printf(flags->mensaje);
        }
         if (strncmp(buffer, "Mostrar tiempos", 15) == 0) {
                flags->mostrar_tiempos = 1;
                flags->servidor_ocupado = 1;
                while(flags->servidor_ocupado==1)
                {
					printf("Esperando respuesta ...\n");
					sleep(3);
				}
				for(id=0;flags->tiempos[id][1]!=-1;id++)
				{
					printf("ID: %d, en T=%d\n",id,flags->tiempos[id][0]);
				}
        }
        if(strncmp(buffer, "Salir del loop principal",24) == 0){
			flags->salida_loop=1;
			flags->servidor_ocupado=1;
		}
		if(strncmp(buffer, "Matar al proceso",16) == 0){
			flags->salida_violenta=1;
			flags->servidor_ocupado=1;
		}
		if (strncmp(buffer, "Crono:", 6) == 0) {
			sscanf(buffer,"Crono:%d",&Tiempo);
			alarm(Tiempo);
			FLcrono_establecido = 1;
			printf("Crono establecido: Guarda estado cada %d segundos\n",Tiempo);     
        }
         if (strncmp(buffer, "Elimina Crono", 13) == 0) {
               alarm(0);
        }
		if(strncmp(buffer,"Reconectarse",12)==0)
		{
			if (shmdt(mem) == -1) {
			fprintf(stderr, "shmdt fallo\n");  
			}else{
				printf("Clave de conexion:\n");
				fgets(buffer, 100, stdin);
				sscanf(buffer, "%d", &clave);
				printf("%d\n",clave);

				if((shmid = shmget((key_t)clave, sizeof(interfase), 0660 | IPC_CREAT)) == -1)
				{
					fprintf(stderr, "shmget fallo\n");
				}

				if((mem = shmat(shmid, NULL, 0))==(void *)-1)
				{
					fprintf(stderr, "shmat fallo\n");
				}else{
					printf("Enlazado al proceso con clave:%d\n",clave);
				
					flags = (interfase *)mem;
				}
			}
		}
        if (strncmp(buffer, "salir", 5) == 0) {
                running = 0;
        }
    }

    if (shmdt(mem) == -1) {
        fprintf(stderr, "shmdt fallo\n");
        exit(EXIT_FAILURE);
    }
    exit(EXIT_SUCCESS);
}

void GuardaEstado(void)
{
	
return;
}

