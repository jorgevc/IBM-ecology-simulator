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
#include <sys/stat.h>
#include <stdlib.h>
//#include "libPP_3.0.h"
//#include "libPP_4.0.h"
#include "libPP_5.0.h"
#include "EntSalArb_MP.h"

void GuardaEstado(estado *es, FILE *archivo)
{
int **s=es->s;
int **TIPO=es->TIPO;
int NDX = es->NDX;
int NDY = es->NDY;
int i,j;


	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]!=0){ fprintf(archivo,"%d   %d   %d\n",i,j,TIPO[i][j]); }
			//else{fprintf(archivo,"No hay datos para: %d   %d \n",i,j);}
		}
	}
}

FILE* OpenFile(char *nombre)
{
FILE *aA;
char archivo[200]="DATOS/";
// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.
strcat(archivo,nombre);
		aA=fopen(archivo, "w");
		fputs("# x   y   tipo\n",aA);
		return aA;
}

void CreaContenedor(char *nombre)
{
char dir[200]="DATOS/";

strcat(dir,nombre);
mkdir(dir,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH));

fprintf(stdout,"Contenedor creado:\n %s \n",nombre);
return;
}

void GuardaEstadoEn(char *nombre, estado *es)
{
int T=es->T;
char paso[15];
char archivo[200];
FILE *datos;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;

sprintf(paso,"/T_%03d",T);
strcpy(archivo,nombre);
strcat(archivo,paso);

	datos=OpenFile(archivo);
	GuardaEstado(es, datos);
	fclose(datos);
	
return;
}


FILE* AbreRhoVsTEn(char *contenedor)
{
FILE *aA;
char archivo[200]="DATOS/";
// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.
strcat(archivo,contenedor);
strcat(archivo,"/RhoVsT");
		aA=fopen(archivo, "w");
		fputs("# t   rho   tipo (tipo 0 es la total)\n",aA);
		fclose(aA);
		aA=fopen(archivo, "a");
		return aA;
}

float ActualizaRhoVsT(estado *es,FILE *archivo,int NoEspecies)
{
int T=es->T;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
int **TIPO=es->TIPO;
sitio *SO=es->SO;
float rho,rho_specie;
int tot=(NoEspecies+1);
int rhoVec[tot];  
memset(rhoVec,0,tot * sizeof(int));

rho=((float)ON)/((float)(NDX*NDY));

fprintf(archivo,"%d   %f   0\n",T,rho); // Tipo 0 es el total

	if(NoEspecies!=0)
	{
	int n;
		for(n=1;n<=ON;n++)
		{
			rhoVec[TIPO[SO[n].i][SO[n].j]]+=1;
		}

		for(n=1;n<tot;n++)
		{
				if(rhoVec[n]!=0)
				{
					rho_specie=((float)rhoVec[n])/((float)(NDX*NDY));
					fprintf(archivo,"%d   %f   %d\n",T,rho_specie,n);
				}
		}
	}
	
return;
}

int GuardaTiposEn(char *contenedor, estado *es)
{
int T=es->T;
char paso[25];
char archivo[200]="DATOS/";
FILE *datos;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
int **s=es->s;
int **TIPO=es->TIPO;
int i,j;
int tot=NDX*NDY;
int rho[tot];  // CUIDADO!: NO PUEDE HABER ETIQUETAS DE ESPECIES MAS GRANDES QUE EL NUMERO DE SITIOS EN LA RED !!!! SEGMENTATION FAULT!!!! 
memset(rho,0,tot * sizeof(int));

sprintf(paso,"/SpeciesRhoRankT_%03d",T);
strcat(archivo,contenedor);
strcat(archivo,paso);

	// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.
	datos=fopen(archivo, "w");
	fputs("# rank   Rho\n",datos);

		for(i=1;i<=NDX;i++)
		{
			for(j=1;j<=NDY;j++)
			{
				if(s[i][j]!=0){ 
				rho[TIPO[i][j]]+=1;
				 }
			}
		}
		
		int tmp;
		int TOTAL=tot-1;
		int *rank=rho;
		i=0;
		int NoSpecies=0;
		
			while(i<TOTAL)
			{
				if(rank[i]<rank[i+1])
				{
					tmp=rank[i];
					rank[i]=rank[i+1];
					rank[i+1]=tmp;
					if(i>0){i--;}
				}else{
				i++;
				}	
			}
		
	float rang_relativo;
		for(i=1;i<tot;i++)
		{
			if(rank[i]!=0)
			{
				rang_relativo=((float)rank[i])/((float)ON);
				fprintf(datos,"%d   %f\n",i,rang_relativo);
				NoSpecies+=1;
			}
		}
	
	fclose(datos);
	
return NoSpecies;
}

FILE* AbreNoSpeciesVsTEn(char *contenedor)
{
FILE *aA;
char archivo[200]="DATOS/";
// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.
strcat(archivo,contenedor);
strcat(archivo,"/SpeciesVsT");
		aA=fopen(archivo, "w");
		fputs("# t   No_Species \n",aA);
		fclose(aA);
		aA=fopen(archivo, "a");
		return aA;
}

void ActualizaNoSpeciesVsT(FILE *archivo,int Species, int T)
{
fprintf(archivo,"%d   %d\n",T,Species);
return;
}


void GuardaCorrelacion(estado *es,int Rini, int Rfin,char *contenedor)
{
FILE *corr;
int T=es->T;
int r;
float g;
char archivo[250]="DATOS/";
char nombre[200];
sprintf(nombre,"/Correlacion_%03d",T);
strcat(archivo,contenedor);
strcat(archivo,nombre);

        corr=fopen(archivo,"w");
	if(corr==NULL){
		puts("No se pudo abrir archivo para guardar correlacion\n");
		return;
	}

	fputs("#r   g\n",corr);
	for(r=Rini;r<=Rfin;r++)
	{
		g=Correlacion(es,r);
		fprintf(corr,"%d   %f\n",r,g);
	}
        fclose(corr);
        
return;
}

void GuardaCorrelacionTipo(estado *es,int Rini, int Rfin,int TOrigen,int TObjetivo,char *contenedor)
{
FILE *corr;
int T=es->T;
int r;
float g;
char archivo[250]="DATOS/";
char nombre[200];

sprintf(nombre,"/CorrelacionT_%03d",T);
strcat(archivo,contenedor);
strcat(archivo,nombre);
puts(archivo);
        corr=fopen(archivo,"w");
	if(corr==NULL){puts("No se pudo abrir archivo");}
	fputs("#r   g\n",corr);
	printf("calculando correlacion al tiempo T=%d de r_ini=%d r_fin=%d\n",T,Rini, Rfin);
	for(r=Rini;r<=Rfin;r++)
	{
		g=CorrelacionEspecies(es,r,TOrigen, TObjetivo);
		fprintf(corr,"%d   %f\n",r,g);
	}
        fclose(corr);
        
return;
}

void GuardaCorrelacion_MP(char *contenedor, Float1D_MP *corr)
{
	FILE *arch;
	char archivo[250]="DATOS/";
	char nombre[200];
	int Rfin=corr->i_max;
	int T=corr->T;
	int r;
	
	sprintf(nombre,"/CorrelacionT_%03d",T);
	strcat(archivo,contenedor);
	strcat(archivo,nombre);
	
	arch=fopen(archivo,"w");
	if(arch==NULL){puts("No se pudo abrir archivo");}
	fputs("#r   g\n",arch);
	for(r=1;r<=Rfin;r++)
	{
		fprintf(arch,"%d %f\n",r,corr->array[r]/corr->NoEnsambles);
	}
        fclose(arch);
        
return;
}

void GuardaCorrelacionTipo_MP(char *contenedor, Float1D_MP *corr)
{
	FILE *arch;
	char archivo[250]="DATOS/";
	char nombre[200];
	int Rfin=corr->i_max;
	int T=corr->T;
	int r;
	
	sprintf(nombre,"/CorrelacionTipoT_%03d",T);
	strcat(archivo,contenedor);
	strcat(archivo,nombre);
	
	arch=fopen(archivo,"w");
	if(arch==NULL){puts("No se pudo abrir archivo");}
	fputs("#r   g\n",arch);
	for(r=1;r<=Rfin;r++)
	{
		fprintf(arch,"%d %f\n",r,corr->array[r]/corr->NoEnsambles);
	}
        fclose(arch);
        
return;
}

int CargaEstado(char *contenedor, char *nombre, estado *es,int NDX, int NDY)   //Usar sin haber alojado memoria antes!!!
{
FILE *datos=NULL;
int **s;
int **TIPO;
sitio *SO;
int **INDICE;
int Tiempo;
int i,j,t;
int max_i=1;
int max_j=1;
int n=0;
char *buffer;
size_t tam_buffer=100*sizeof(char);
 buffer = (char *) malloc (tam_buffer + 1);
int args_assigned = 0;

char archivo[250]="DATOS/";

strcat(archivo,contenedor);
strcat(archivo,"/");
strcat(archivo,nombre);

puts(archivo);

	if((datos = fopen (archivo, "r"))==NULL){
		puts("\nNo se pudo abrir para leer\n");
		return 0;
		}
	
	while(getline(&buffer, &tam_buffer, datos)!=-1)
	{
		if(strchr(buffer, '#')==NULL)
		{
			args_assigned = sscanf(buffer, "%d %d", &i, &j);
			if(i>max_i){max_i=i;}
			if(j>max_j){max_j=j;}
		}
	}
	
	if(max_i<NDX){max_i=NDX;}
	if(max_j<NDY){max_j=NDY;}
	
	AlojaMemoria(es,max_i,max_j);
	ResetEstado(es);
	s=es->s;
	TIPO=es->TIPO;
	SO=es->SO;
	INDICE=es->INDICE;
	
	rewind( datos );
	puts("leyendo...\n");
		
	while (getline(&buffer, &tam_buffer, datos)!=-1)
    {
      args_assigned = sscanf (buffer, "%d %d %d", &i, &j, &t);
      if(args_assigned == 3)
      {
		  n++;
		s[i][j]=1;
		TIPO[i][j]=t;
		SO[n].i=i;
		SO[n].j=j;
		INDICE[i][j]=n;
		(es->ON)++;  
	  }else{
		  if(args_assigned == 2)
		  {
			   n++;
			s[i][j]=1;
			TIPO[i][j]=0;
			SO[n].i=i;
			SO[n].j=j;
			INDICE[i][j]=n;
			(es->ON)++;  
		}
	  }  
    }
    
    fclose(datos);
    printf("leidas %d lineas\n",n);
    
    if(sscanf(nombre,"T_%d",&Tiempo)==1)
    {
		es->T=Tiempo;
	}else{
		es->T=0;
	}
	
	printf("Tiempo asignado: T=%d\n",es->T);
	
return 1;
	
}

void GuardaRhoVsT_MP(char *contenedor, Float2D_MP *RhoVsT, Dist_MP *RhoDist)
{
FILE *datos, *dist;

if(RhoVsT!=NULL)
{
	int T_max=RhoVsT->i_max;
	int NoEspecies=RhoVsT->j_max;
	float NoEnsambles=(float)RhoVsT->NoEnsambles;

	datos=AbreRhoVsTEn(contenedor); 

	int T,e;
	for(T=0;T<=T_max;T++)
	{
		for(e=0;e<=NoEspecies;e++)
		{
				fprintf(datos,"%d %f %d\n",T,RhoVsT->array[T][e]/NoEnsambles,e);
		}
		
	}
	fclose(datos);
}

if(RhoDist!=NULL)
{
	int RhoPart;
	char archDist[250]="DATOS/";
	char nombre[30];
	strcat(archDist,contenedor);
	sprintf(nombre,"/RhoDistT_%d",RhoDist->T);
	strcat(archDist,nombre);
	dist=fopen(archDist,"w");
	fputs("# Rho   Prob   T\n",dist);
		for(RhoPart=0;RhoPart<=(RhoDist->i_max);RhoPart++)
		{
			fprintf(dist,"%f %f\n",((float)RhoPart) * (RhoDist->TamParticion),(float)(RhoDist->array[RhoPart])/(float)(RhoDist->NoEnsambles));
		}
		
	fclose(dist);
}

return;
}

void GuardaTiposEn_MP(char *contenedor,Float2D_MP *MP_RhoVsT,int T)
{	
	char paso[25];
	char archivo[200]="DATOS/";
		FILE *datos;
	    int i;
	    int MaxEspecie=MP_RhoVsT->j_max;
		float rank[MaxEspecie + 1];
		int NoSpecies=0;
		float NoEnsambles=(float)MP_RhoVsT->NoEnsambles;
		float tmp;
		
		for(i=0;i<=MaxEspecie;i++)
		{
			rank[i]=MP_RhoVsT->array[T][i];
		}
		
			i=0;
			while(i<MaxEspecie)
			{
				if(rank[i]<rank[i+1])
				{
					tmp=rank[i];
					rank[i]=rank[i+1];
					rank[i+1]=tmp;
					if(i>0){i--;}
				}else{
				i++;
				}	
			}
			
		sprintf(paso,"/SpeciesRhoRankT_%03d",T);
		strcat(archivo,contenedor);
		strcat(archivo,paso);
		
		datos=fopen(archivo,"w");
		float rang_relativo;
		for(i=1;i<=MaxEspecie;i++)
		{
			if(rank[i]!=0)
			{
				rang_relativo=(rank[i])/(NoEnsambles);
				fprintf(datos,"%d   %f\n",i,rang_relativo);
				NoSpecies+=1;
			}
		}
		fclose(datos);
		
return;
}

void GuardaEstadoEn_MP(char *nombre, estado *es,int id,int ensamble)
{
int T=es->T;
char paso[15];
char archivo[300];
char base[50];
FILE *datos;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;

char dir[250]="DATOS/";

sprintf(paso,"/T_%03d",T);
strcat(dir,nombre);
strcat(dir,paso);
mkdir(dir,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH));

sprintf(base,"/P_%d_Ens_%d",id,ensamble);
strcpy(archivo,dir);
strcat(archivo,base);

	datos=fopen(archivo,"w");
	fputs("# x   y   tipo\n",datos);
	GuardaEstado(es, datos);
	fclose(datos);
	
return;
}

int CargaEstado_MP(char *contenedor, char *nombre, estado *es,int NDX, int NDY,int id, int NoEnsambles)  //Usar sin haber alojado memoria antes!!!
{

char ensamble[25];
char contenedor_MP[250];
int ens;
int Tiempo;
int cargados=0;
strcpy(contenedor_MP,contenedor);
strcat(contenedor_MP,"/");
strcat(contenedor_MP,nombre);

		if(sscanf(nombre,"T_%d",&Tiempo)!=1)
		{
			Tiempo=0;
		}

int status=1;
	for(ens=0;status==1;ens++)
	{
		sprintf(ensamble,"P_%d_Ens_%d",id,ens);
		status=CargaEstado(contenedor_MP,ensamble,&es[ens],NDX,NDY);	
		if(status==1)
		{
			es[ens].T=Tiempo;
			printf("Tiempo reasignado %d\n",Tiempo);
			cargados++;
		}
		if(cargados>=NoEnsambles)
		{
			status=-1;
		}
	}

return cargados;	
}
