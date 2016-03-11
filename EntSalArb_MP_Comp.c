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
#include "libPP_6.1.h"
#include "EntSalArb_MP_Comp.h"
#include <errno.h>

void GuardaEstado(estado *es, FILE *archivo)
{
int **s=es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i,j;


	for(i=1;i<=NDX;i++)
	{
		for(j=1;j<=NDY;j++)
		{
			if(s[i][j]>0){ 
				fprintf(archivo,"%d %d %d %d\n",i,j,es->individuals[es->INDICE[i][j]].species,(int)es->individuals[es->INDICE[i][j]].size_float); 
				}
			//else{fprintf(archivo,"No hay datos para: %d   %d \n",i,j);}
		}
	}
}


void CreaContenedor(char *nombre,runDescriptor run)
{
	fprintf(stdout,"\nCreando contenedor: %s",nombre);
	
char copyFileName[550];
char dir[550];
char *pch;
FILE *file;

strcpy(dir,"");
	strcpy(copyFileName,nombre);
	pch = strtok (copyFileName,"/");
	while(pch != NULL)
	{	
		strcat (dir,pch);
		 strcat (dir,"/");
		if(mkdir(dir,(S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) != 0)
		 {
			 fprintf(stdout,"\ndir: %s, ",dir);
			 perror("Error creando directorio: ");
		 }
		 pch = strtok (NULL, "/");
	}

strcpy(copyFileName,nombre);
strcat(copyFileName,"/descriptor");
file=fopen(copyFileName,"w");

#ifdef EXPLICIT_RESOURCES
fprintf(file,"Algorithm: EXPLICIT_RESOURCES\n\n");
#endif
#ifdef VIRTUAL_GRID
fprintf(file,"Algorithm: VIRTUAL RESOURCES GRID\n\n");
#endif
#ifdef SOI
fprintf(file,"Algorithm: SOI\n\n");
#endif

fprintf(file,"X=%d\n",run.X);
fprintf(file,"Y=%d\n",run.Y);
fprintf(file,"grid_units=%f  (Area_units=units^2)\n",run.grid_units);
fprintf(file,"size_units=%f",run.size_units);
fprintf(file,"T_max=%d\n",run.T_max);
fprintf(file,"NoEnsambles=%d\n", run.NoEnsambles);

fprintf(file,"\ncoagulation_exp=%f (Proportional to the Area)\n", run.Model.coagulation_exp);
fprintf(file,"coagulation_factor=%f  (Area_units/(Area_units)^c_exp)*factor\n", run.Model.coagulation_factor);

fprintf(file,"\ncoagulation_radio_exp=%f\n", run.Model.coagulation_radio_exp);
fprintf(file,"coagulation_radio_factor=%f  (units/(Area_units)^cr_exp)*factor\n", run.Model.coagulation_radio_factor);

fprintf(file,"\ncoagulation_metabolic_exp=%f\n", run.Model.metabolic_exp);
fprintf(file,"metabolic_factor=%f  (Area_units/(Area_units)^m_exp)*factor\n", run.Model.metabolic_factor);

fprintf(file,"\ngrowth_constant=%d  ([int]>0 needed resources per unit size increse.)",run.Model.growth_constant);

fprintf(file,"\nhealth_factor=%f\n", run.Model.health_factor);
fprintf(file,"\nmin_health=%d\n", run.Model.min_health);

fprintf(file,"\nresource_rate=%f\n", run.Model.resource_rate);
fprintf(file,"ResourcesScale=%d\n", run.Model.ResourcesScale);
fprintf(file,"competitionAsymetry=%f\n", run.Model.competitionAsymetry);

fprintf(file,"\nbirth_rate=%f\n",run.Model.birth_rate);
fprintf(file,"RadioBirth=%d\n",run.Model.RadioBirth);
fprintf(file,"dead_rate=%f (intrinsic dead rate)\n", run.Model.dead_rate);

fprintf(file,"---INITIAL CONDITIONS---\n");
fprintf(file,"Initial mean distance among individuals=%d\n",run.initialMeanDistance);
fprintf(file,"Initial minimun separation among individuals=%d\n", run.initialMinSeparation);

fclose(file);
	
fprintf(stdout,"\nlisto creando contenedor!\n");


return;
}

void GuardaEstadoEn(char *contenedor, estado *es)
{
int T=es->T;
char archivo[550];
FILE *datos;

sprintf(archivo,"%s/T_%03d",contenedor,T);

	datos=fopen(archivo, "w");
	fputs("# x y tipo tamano\n",datos);
	GuardaEstado(es, datos);
	fclose(datos);
	
return;
}


FILE* AbreRhoVsTEn(char *contenedor)
{
FILE *aA;
char archivo[200];
// mkdir(archivo,00777); Descomentar para crear el directorio si es que no existe.

sprintf(archivo,"%s/RhoVsT",contenedor);
		aA=fopen(archivo, "w");
		fputs("# t   rho   tipo (tipo 0 es la total)\n",aA);
		return aA;
}

float ActualizaRhoVsT(estado *es,FILE *archivo,int NoEspecies)  //Deprecated
{
int T=es->T;
int NDX = es->NDX;
int NDY = es->NDY;
int ON = es->ON;
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
			rhoVec[es->individuals[n].species]+=1;
			//rhoVec[s[SO[n].i][SO[n].j]]+=1;
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
				rho[es->individuals[es->INDICE[i][j]].species]+=1;
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


void GuardaCorrelacion_MP(char *contenedor, char *prefix, Float1D_MP *corr)
{
	FILE *arch;
	char archivo[250]="DATOS/";
	char nombre[200];
	int Rfin=corr->i_max;
	
	int T=corr->T;
	int r;
	
	sprintf(nombre,"/%s_CorrT_%03d",prefix,T);
	strcat(archivo,contenedor);
	strcat(archivo,nombre);
	
	arch=fopen(archivo,"w");
	if(arch==NULL){puts("No se pudo abrir archivo");}
	fputs("#r g\n",arch);
	for(r=1;r<=Rfin;r++)
	{
		//if(corr->array[r]!=0.0)
		//{
			fprintf(arch,"%d %f\n",r,corr->array[r]/((float)corr->NoEnsambles));
		//}
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
	fputs("#r g\n",arch);
	for(r=1;r<=Rfin;r++)
	{
		if(corr->array[r]!=0.0)
		{
			fprintf(arch,"%d %f\n",r,corr->array[r]/((float)corr->NoEnsambles));
		}
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
	
	SO=es->SO;
	INDICE=es->INDICE;
	
	rewind( datos );
	puts("leyendo...\n");
		
	while (getline(&buffer, &tam_buffer, datos)!=-1)
    {
      args_assigned = sscanf (buffer, "%d %d %d", &i, &j, &t);
      if(args_assigned == 3)
      { 
		InsertaIndividuoEn(es,i,j,t,1);
	  }else{
		  if(args_assigned == 2)
		  {
			 InsertaIndividuoEn(es,i,j,0,1);
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
				//printf("guadando:%d %f %d NoEnsambles=%d \n",T,RhoVsT->array[T][e]/NoEnsambles,e, RhoVsT->NoEnsambles);
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

void GuardaCorrXY(Float2D_MP *correlacion, char *contenedor,char *sufix)
{
	FILE *Arch;
	int i,j;
	int NDX = correlacion->i_max;
	int NDY = correlacion->j_max;
	char nombre[150];
	
		printf("Guardando CFFT_XY\n");
		sprintf(nombre,"DATOS/%s/CFFT_%s_MP_XY",contenedor,sufix);
		Arch = fopen(nombre,"w");
		for(i=0;i<NDX;i++)
		{
			for(j=0;j<NDY;j++)
			{
				if(correlacion->array[i][j]>0.0)
				{
					fprintf(Arch,"%d %d %f\n",i,j,correlacion->array[i][j]/((float)correlacion->NoEnsambles));
				}
			}
		}
		fclose(Arch);
		printf("Se ha Guardado CFFT_XY\n");
return;
}

void GuardaFloat1D_MP(char *contenedor,char *nombre, Float1D_MP *MP_Float1D)
{
	FILE *Arch;
	int i;
	int NDX = MP_Float1D->i_max;
	float delta_s=1.0/MP_Float1D->index_units;
	char nombreCom[150];
	
		printf("Guardando .. ");
		sprintf(nombreCom,"%s/%s",contenedor,nombre);
		printf("\nen %s",nombreCom);
		Arch = fopen(nombreCom,"w");
		for(i=0;i<=NDX;i++)
		{
				if(MP_Float1D->array[i]!=0.0)
				{
					fprintf(Arch,"%f %f\n",((float)i)*delta_s,MP_Float1D->array[i]/((float)MP_Float1D->NoEnsambles));
				}
		}
		fclose(Arch);
		printf("Se ha Guardado.\n");
return;
	
}

void GuardaDist_MP(char *contenedor,char *nombre, Dist_MP *MP_Dist)
{
	FILE *dist;
int RhoPart;
	char archDist[250]="DATOS/";
	strcat(archDist,contenedor);
	strcat(archDist,"/");
	strcat(archDist,nombre);
	dist=fopen(archDist,"w");
	fputs("# x Prob\n",dist);
		for(RhoPart=0;RhoPart<=(MP_Dist->i_max);RhoPart++)
		{
			if(MP_Dist->array[RhoPart]!=0)
			{
				fprintf(dist,"%f %f\n", (MP_Dist->xIni + ((float)RhoPart) * (MP_Dist->TamParticion)),(float)(MP_Dist->array[RhoPart])/(float)(MP_Dist->NoEnsambles));
			}
		}
		
	fclose(dist);
	
	return;
}
