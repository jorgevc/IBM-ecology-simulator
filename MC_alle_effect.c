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

#include "MC_alle_effect.h"
#include "GNA.h"
#include <stdio.h>
#include <time.h>
#include <sqlite3.h>

void MC_sweep_alle(estado *es, alle_env *env)
{
//	start2 = clock();  //clock comentar!!
int Indice;
float DT=0.0;
int Total_Individuals=es->ON;
	
	while(DT<1.0){
		if(Total_Individuals != 0){
			DT+=1.0/Total_Individuals; 
			Indice = I_JKISS(1,Total_Individuals);
			Update_alle(es, env, Indice);
			Total_Individuals = es->ON;	
		}else{
			DT=2.0;
		}
	}
	
(es->T)++;		

//     end2 = clock();			//clock comentar!!
  //   tiempo2 = ((double) (end2 - start2))/ CLOCKS_PER_SEC;;  //clock comentar!!	
  //   printf("Tiempo en elejir vecinos en una corrida:%f, tiempo corrida:%f\n",cpu_time_used,tiempo2);
   //  cpu_time_used=0.0;
		
return;
}


void Update_alle(estado *es, alle_env *env, int N)
{
float Max_Metabolic = env->Max_Metabolic;
float Rand; 
float pDead, pCreacion, pCoagulation1, C, Dead, Birth, Coagulation1, Nada,pCoa2;
sitio vecino;
int **s = es->s;
int NDX = es->NDX;
int NDY = es->NDY;
int i=es->SO[N].i;
int j=es->SO[N].j;
sitio *SO = es->SO;
int radioCre;
int radioCoa;
int ho = (es->TIPO[i][j] * 10) - 5;	//Nicho0 cambiar cuando cambie numero de especies distinto a 50 (solo en nicho es necesario)
	
	Rand = F_JKISS();
	
	Dead=env->param[es->TIPO[i][j]].Dead; // modelo neutral y J-C
	// rdead = ((float)(es->TIPO[i][j]-500)/5000.0)+0.05; //Modelo Lottery
	Birth=env->param[es->TIPO[i][j]].Birth;
//	Dead=Birth - (Birth - parametros[es->TIPO[i][j]].Dead) * exp(-pow(campo - ho,2.0)/20000.0);  // funcion nicho0 remaster
	Coagulation1=env->param[es->TIPO[i][j]].Coagulation;	
	pCoa2=env->param[es->TIPO[i][j]].CoagulationIntra/Max_Metabolic;
	pDead=Dead/Max_Metabolic;	//Asignar Max_Metabolic, si no hay division entre cero.
	pCreacion = Birth/Max_Metabolic;
	pCoagulation1 = Coagulation1/Max_Metabolic;

	if(Rand<=pDead) //aniquilacion
	{	
		s[i][j]=0;
		SO[N]=SO[(es->ON)];
		es->INDICE[SO[es->ON].i][SO[es->ON].j]=N;
		(es->ON)--;	
		es->TIPO[i][j]=0;
		
	}else{  //creation o coagulacion o nada
		if(Rand<=(pDead + pCreacion)) //creation
		{
			int alle_go=1;
			if(env->alle_effect > 0)
			{
				radioCre=env->param[es->TIPO[i][j]].alle_range;
				EligeUniforme(i,j,radioCre,&vecino);
				
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;} 
				if(s[vecino.i][vecino.j] <= 0)
				{
					if(Rand<=(pDead + pCreacion*env->alle_effect))
					{
						alle_go=0;
					}
				}
			}
			
			if(alle_go==1)
			{
				radioCre=env->param[es->TIPO[i][j]].RadioBirth;
				EligeUniforme(i,j,radioCre,&vecino);
				
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
				
				if(s[vecino.i][vecino.j]==0)
				{
					s[vecino.i][vecino.j]=1;
					(es->ON)++;
					SO[(es->ON)]=vecino;
					es->INDICE[vecino.i][vecino.j]=(es->ON);
					es->TIPO[vecino.i][vecino.j]=es->TIPO[i][j];
				}
			}
			
		}else{ //coagulacion1 o 2 o nada
			if(Rand<=(pDead + pCreacion + pCoagulation1))  //coagulacion
			{
				radioCoa=env->param[es->TIPO[i][j]].RadioCoa;
				 EligeUniforme(i,j,radioCoa,&vecino);
				if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
				if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
				if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
				if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
				
				if(s[vecino.i][vecino.j]==1 && (es->TIPO[vecino.i][vecino.j])!=(es->TIPO[i][j]))     // mato a los que NO SON de mi propia especie 
				{
					s[vecino.i][vecino.j]=0;
					SO[(es->INDICE[vecino.i][vecino.j])]=SO[(es->ON)];
					es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[vecino.i][vecino.j]);
					(es->ON)--;
					es->TIPO[vecino.i][vecino.j]=0;
				}
			}else{  //coagulacion 2 o nada
				if(Rand<=(pDead + pCreacion + pCoagulation1 + pCoa2)) //coagulation Inter
				{
					radioCoa=env->param[es->TIPO[i][j]].RadioCoaIntra;
					 EligeUniforme(i,j,radioCoa,&vecino);
					if(vecino.i <= 0){vecino.i = NDX + vecino.i;}
					if(vecino.j <= 0){vecino.j = NDY + vecino.j;}
					if(vecino.i > NDX){vecino.i = vecino.i - NDX;}
					if(vecino.j > NDY){vecino.j = vecino.j - NDY;}   //NOTA: Peligro de segmentation fault si el radio es mayor al lado de la maya
					
					if(s[vecino.i][vecino.j]==1 && (es->TIPO[vecino.i][vecino.j])==(es->TIPO[i][j]))     //  Solo mato a mi propia especie 
					{
						s[vecino.i][vecino.j]=0;
						SO[(es->INDICE[vecino.i][vecino.j])]=SO[(es->ON)];
						es->INDICE[SO[es->ON].i][SO[es->ON].j]=(es->INDICE[vecino.i][vecino.j]);
						(es->ON)--;
						es->TIPO[vecino.i][vecino.j]=0;
					}
				}	
			}
		}
	}
return;
}

float calculate_metabolic_time(alle_env *env)
{	
	float MaxMetabolic=0.0;
	if(env->param==NULL)
	{
		printf("No hay parametros a normalizar\n");
		return;
	}else{
		int i;
		float Metabolizmo=0.0;
		
		for(i=1;i<=env->NumberOfSpecies;i++)
		{
			Metabolizmo=env->param[i].Birth;
			Metabolizmo+=env->param[i].Coagulation;
			Metabolizmo+=env->param[i].CoagulationIntra;
			Metabolizmo+=env->param[i].Dead;
			if(MaxMetabolic<Metabolizmo)
			{
				MaxMetabolic=Metabolizmo;
			}
		}
	}
return MaxMetabolic;
}

int store_sim_db(char *values)
{
	sqlite3 *db;
    char *err_msg = 0;
    int inserted_id = -1;
    int rc = sqlite3_open("alle.db", &db);
    if( rc ){
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
   }else{
		char sql[5000];
		sprintf(sql, "INSERT INTO single_species VALUES (%s);", values);
		rc = sqlite3_exec(db, sql, 0, 0, &err_msg);
		inserted_id = sqlite3_last_insert_rowid(db);
		sqlite3_close(db);
	}
return inserted_id;
}
