/*
Copyright 2012 Jorge Velazquez
*/
typedef struct {
	int noPoints;
	int *point74;
	int *point84;
	int *point93;
	int noPoints84;
	int noPoints93;
} DataSet;

typedef struct {
	int A;
	int B;
} Pair;


Pair meanDeviation(int *data, int noPoints, estado *e, float scale)
{
	int state[noPoints];
	int total=0;
	int i;
	int distance=0;
	Pair deviation;
	deviation.A=0;
	deviation.B=0;
	
	for(i=0;i<noPoints;i++)
		{
			state[i]=0;
		}
			
		for(i=1;i<=e->ON;i++)
		{
			state[e->individuals[i].index]=scale*(e->individuals[i].size_float);
		}
		
		for(i=0;i<noPoints;i++)
		{
			if(data[i]>0)
			{
				deviation.A -= (state[i] - data[i])*(state[i] - data[i]);
				deviation.B++;
			}
		}
		
		return deviation;	
}

Pair HowIsDistSizesSymilar(int *data,int noPoints,estado *e, int delta, float scale)
{
	int state[noPoints];
	int ON=e->ON;
	int i;
	Pair like;
		
		for(i=0;i<noPoints;i++)
		{
			state[i]=0;
		}
			
		for(i=1;i<=ON;i++)
		{
			state[e->individuals[i].index]=scale*(e->individuals[i].size_float);
		}
		
	int	coincidences=0;
			for(i=0;i<noPoints;i++)
			{
					if((data[i] < (state[i] + delta)) && (data[i] > (state[i] - delta)))
					{
						coincidences++;
					}
			}
		
		like.B=noPoints;
		like.A=coincidences;
		
	return like;
}

Pair MLE(DataSet *data,runDescriptor run, int *time_step84, int *time_step93)
{	

int const write_interval=(run.grid_units*run.grid_units)*200;

////
int T_max = run.T_max;
int NoEnsambles=run.NoEnsambles;


model modelo;
modelo = run.Model;

 

int NDX=(run.X*run.grid_units)+1;
int NDY=(run.Y*run.grid_units)+1;

omp_set_num_threads(4);

////////////////////////////Termina Inicializa de paremtros de la simulacion
/////////////////////////////////////Prepara CONTENEDOR para escribir DATOS:
//Float1D_MP TamDist;
//	InicializaFloat1D_MP(&TamDist, T_max+100);
//char distT[10];
//char contenedor[150];
//	sprintf(contenedor,"DATOS_MCMC/prueba");


Pair GTotalmetaSymilar84;
Pair GTotalmetaSymilar93;
int timeOfLike84;
int timeOfLike93;

GTotalmetaSymilar84.A=0;
GTotalmetaSymilar84.B=0;

float years=10.0;
///////////////////////////////////// Estado INICIAL:

	#pragma omp parallel			///////INICIA PARALLEL
	{	
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
			e[Par].units=run.grid_units;
			e[Par].size_units=run.size_units;
		}

		init_JKISS(); //Inicializa la semilla de cada proceso.
float scale=run.scaleFactor;	


Individual indv;
indv.species=1;
indv.size_float=5;
indv.size=run.size_units*10.0*indv.size_float;
indv.radio_float=R(indv, (&modelo));
indv.radio=indv.radio_float;
indv.metabolism=0;
indv.health=0;

int flag[MaxPar+1];
flag[MaxPar]=0;
int i,j,n;
			for(Par=0;Par<MaxPar;Par++)
			{
				while(e[Par].ON < data->noPoints)
				{
					i=I_JKISS(1,run.X);
					j=I_JKISS(1,run.Y);
					indv.size_float=((float)data->point74[e[Par].ON])/scale;
					indv.radio_float=R(indv, (&modelo));
					indv.index=e[Par].ON;
					InsertIndividualAt(&e[Par],i,j,indv,0);				
				}
			setMaxMetabolic(&e[Par],&modelo);	
			flag[Par]=0;
			}		
	/////////////////////////////////////Termina Estado INICIAL
	//////////////////////////Prepara Contenedor en Memoria de cada proceso Para Mejorar rendimiento (optimizar el uso de cache de cada procesador)

	//////////////////////////////Barrido Monte CARLO:


Pair similar,TotalmetaSymilar84,TotalmetaSymilar93;

TotalmetaSymilar84.A=0;
TotalmetaSymilar84.B=0;

		for(i=1;i<=T_max;i++)
		{			
			for(Par=0;Par<MaxPar;Par++)
			{
				BarrMCcRyCampTamanoSimple(&e[Par], run.Model.resource_rate, &modelo);
			//	ActualizaDistTamano_MP(&e[Par], &TamDist, 'A');
				if(years < e[Par].Meta_T && flag[Par]==0)
				{
					//similar=HowIsDistSizesSymilar(data->point93,data->noPoints,&e[Par],5,scale);
					similar=meanDeviation(data->point93, data->noPoints, &e[Par], scale);
					TotalmetaSymilar84.A+=similar.A;
					TotalmetaSymilar84.B+=similar.B;
					flag[Par]=1;
					flag[MaxPar]++;
					if(flag[MaxPar]==MaxPar)
					{
						i=T_max;		
					}
				}
			}				
			
				if((i-(i/write_interval)*write_interval)==1)    //Inicializa cada write_interval
				{
					init_JKISS();		
				}			
		}

	#pragma omp critical
	{
		GTotalmetaSymilar84.A+=TotalmetaSymilar84.A;
		GTotalmetaSymilar84.B+=TotalmetaSymilar84.B;
	}

		//Libera Memoria
		for(Par=0;Par<MaxPar;Par++)
		{
			LiberaMemoria(&e[Par]);
		}
			
	}	/////TERMINA PARALLEL
	

(*time_step84)=timeOfLike84;
	return GTotalmetaSymilar84;
}

DataSet* loadSizes(char *file, int plot)
{
DataSet *data;
FILE *datos=NULL;
char *buffer;
int dplot,tag,i,next;
size_t tam_buffer = 100*sizeof(char);
buffer = (char *) malloc (tam_buffer);
int args_assigned = 0;
int n=0;
int totalPoints=0;
char pattern[100];

	if((datos = fopen (file, "r"))==NULL){
		puts("\nNo se pudo abrir para leer\n");
		return;
		}
	
	while(getline(&buffer, &tam_buffer, datos)!=-1)
	{
		if(strchr(buffer, '#')==NULL)
		{
			args_assigned = sscanf(buffer, "%d", &dplot );
			if(args_assigned == 1 && plot == dplot)
			{
				totalPoints++;
			}
		}
	}		
	rewind(datos);
	
	data = (DataSet *)malloc(sizeof(DataSet));
	data->point74 = (int *)calloc(totalPoints+1, sizeof(int));
	data->point84 = (int *)calloc(totalPoints+1, sizeof(int));
	data->point93 = (int *)calloc(totalPoints+1, sizeof(int));
	data->noPoints84=0;
	data->noPoints93=0;
	data->noPoints=0;
	
	sprintf(buffer,"000000000000000000000");
	while(getline(&buffer, &tam_buffer, datos)!=-1)
	{
		if(strchr(buffer, '#')==NULL)
		{
			args_assigned = sscanf(buffer, "%d,%d", &dplot, &tag);
			if(args_assigned == 2 && plot == dplot)
			{
				i=0;
					args_assigned = sscanf(buffer, "%*d,%*d,%d", &next);
					if(args_assigned == 1)
					{
						if(next != ',')
						{
							data->point74[n]=next;
							sprintf(pattern,"%%*d,%%*d,%%*d,%%d");
						}else{
							data->point74[n]=0;
							sprintf(pattern,"%%*d,%%*d,,%%d");
							i=1;
						}
						args_assigned = sscanf(buffer, pattern, &next);
						if(args_assigned == 1)
						{
							if(next != ',' && next !=0)
							{
								data->point84[n]=next;
								data->noPoints84++;
								if(i==1)
								{
									sprintf(pattern,"%%*d,%%*d,,%%*d,%%d");
								}else{
									sprintf(pattern,"%%*d,%%*d,%%*d,%%*d,%%d");	
								}
							}else{
								data->point84[n]=0;
								if(i==1)
								{
									sprintf(pattern,"%%*d,%%*d,,,%%d");
								}else{
									sprintf(pattern,"%%*d,%%*d,%%*d,,%%d");
								}	
							}
							args_assigned = sscanf(buffer, pattern, &next);
							if(args_assigned == 1)
							{
								if(next != ',' && next !=0 )
								{
									data->point93[n]=next;
									data->noPoints93++;
								}else{
									data->point93[n]=0;
								}
							}
						}
							
					}
				n++;
			}
		}
		sprintf(buffer,"000000000000000000000");
	}
	
	data->noPoints=n;
	
free(buffer);
fclose(datos);
return data;
}

void FreeDataSet(DataSet *data)
{
	free(data->point74); 
	free(data->point84);
	free(data->point93);
	free(data);
return;
}

int meanDataSet(DataSet *data)
{
	int mean=0;
	int i;
	for(i=0; i < data->noPoints; i++)
	{
		mean+=data->point74[i];
	}
return (mean/data->noPoints);
}
