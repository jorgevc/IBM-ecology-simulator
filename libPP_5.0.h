/** Contains a coordinate pair. 
 * Generaly it is used to label a site on a lattice.
*/
typedef struct {
int i;			/**< i coordinate, could be thinked as the X coordinate*/ 
int j;			/**< j coordinate, could be thinked as the Y coordinate*/ 
} sitio;

/** Contains the parameters describing the individuals of a species. */
typedef struct {
float Birth;			/**< Birth rate of the species*/
float Coagulation;		/**< Inter-Competition rate of the species. 
							Is the rate at wich the precence of an individual of this specie led a dead 
							of an idividual of OTHER species*/
float Dead;				/**< Intrinsic Dead rate.*/
float CoagulationIntra;	/**< Intra-Competition rate.
							Is the rate at wich the precence of an individual of this specie led a dead 
							of an idividual of the SAME species*/
int RadioBirth;			/**< Radio within a descendant could born in units of lattice sites. */
int RadioCoa;			/**< Radio within Inter-Competion acts in units of lattice sites. */
int RadioCoaIntra;		/**< Radio within Intra-Competion acts in units of lattice sites. */
} especie;

/** Contains the state of the system(lattice). */
typedef struct {
int NDX;			/**< The size of the lattice in the i coordinate (X coordinate). */
int NDY;			/**< The size of the lattice in the j coordinate (Y coordinate). */
int T;				/**< The Time-steeps that the sistem has been evolve since the begining of the simulation. 
						It is incremented by one, each time the instance is passed to BarrMCcRyCamp(estado *es) */
int ON;				/**< The total number of occupied sites in the lattice. */
int **s; 			/**< A 2 dimensional array of int representing the actual lattice(system). 
						A value of 0 in s[i][j] represent an empy site in (i,j). A value > 0 represent an occupy site.*/
int **INDICE;		/**< A 2 dimensional array which value at INDICE[i][j] is the index of an OCCUPIED site at (i,j) on a list 
						of the occupied sites. If (i,j) is not occupyed the return value is undetermined. */
sitio *SO;			/**< A 1 dimensional array of scructs sitio that is a list of the occupied sites. 
						To find the index of an specific occupied site at (i,j) in this list, one have to look for it in INDICE[i][j]. */
int **TIPO;			/**< A 2 dimensional array which value at TIPO[i][j] is used to label the species of the individual in (i,j). 
						If the site in (i,j) is not occupyed the return value is undetermined. */
} estado;

/** General porpuse 2 dimensinal array of float suitable to be used for ensemble averages. */
typedef struct {
float **array;		/**< The actual 2 dimenasional array of float: <Float2D_MP>.array[i][j]. */
int i_max;			/**< The maximun value in the i dimension */
int j_max;			/**< The maximun value in the j dimension */
int T;				/**< Used to store the time-step of the property that is stored in the instance(if aplicable). */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
} Float2D_MP;

/** General porpuse 1 dimensinal array of float suitable to be used for ensemble averages. */
typedef struct {
float *array;		/**< The actual 1 dimenasional array of float: <Float1D_MP>.array[i]. */
int i_max;			/**< The maximun value of i. It is equal to (array_size + 1) */
int T;				/**< Used to store the time-step of the property that is stored in the instance(if aplicable). */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
} Float1D_MP;

/** General porpuse 1 dimensinal array of float suitable to be used for ensemble averages. */
typedef struct {
int **array;		/**< The actual 2 dimenasional array of float: <Int2D_MP>.array[i][j]. */
int i_max;			/**< The maximun value in the i dimension */
int j_max;			/**< The maximun value in the j dimension */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
} Int2D_MP;

/** Struct that represent a distribution on a 1 dimensional real space. i.e. Represents a distribution f(x) where x could be real. */
typedef struct {
int *array;			/**< The values of the distribution. */
int T;				/**< Used to store the time-step of the property that is stored in the instance(if aplicable). */
int i_max;			/**< The maximun value of i. It is equal to (array_size + 1) */
float TamParticion; /**< The maximum resolution that the distribution can handle on domain.  */
int NoEnsambles;	/**< The number of ensambles that are stored in the instance. This number is generally used to obtain the mean in a
						multi-thread program, or simplie the ensemble average. */
float xIni;			/**< The initial value of X where the distribution is "defined". */
float xFin;			/**< The final value of X where the distribution is "defined". */
} Dist_MP;

/** Set the value of the Birth rate of a species. 
 * @param[in] L the value of the Birth rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the crated species to save memory. 
 * 			After seting all the values for a specie it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetBirth(float L, int tipo);

/** Set the value of the Inter-Coagulation rate of a species. 
 * @param[in] e the value of the Inter-Coagulation rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero). 
 * 			After seting all the values for a species it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetCoagulation(float e, int tipo);

/** Set the value of the Intra-Coagulation rate of a species. 
 * @param[in] e the value of the Intra-Coagulation rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero). 
 * 			After seting all the values for a species it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetCoagulationIntra(float e,int tipo);

/** Set the value of the Intrinsic Dead rate of a species. 
 * @param[in] d the value of the Intrinsic Dead rate.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero). 
 * 			After seting all the values for a species it should be call the function @see EscalaTiempoMetabolico(int tipo).
 * @see especie
 * */
void SetDead(float d, int tipo);

/** Set the value of the Birth range of a species in number of lattice sites. 
 * @param[in] rb the radio in number of lattice sites of the range of the new offsprings.
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero).
 * @see especie 
 * */
void SetRadioBirth(int rb, int tipo);

/** Set the value of the Inter-Coagulation range of a species in number of lattice sites. 
 * @param[in] rc the radio in number of lattice sites of the range of Inter-Coagulation..
 * @param tipo (>0) the number that identifies a specie. If the specie has not been created, it is created.
 * 			it is recommended to start with values 1,2 ... and so on for the species to save memory (lowest posible numbers greater than zero).
 * @see especie 
 * */
void SetRadioCoa(int rc, int tipo);

/** This shuld be called when a rate of the specie parameters is updated, or after setting up all the rates for the first time.  
 * It sets an internal scale factor of the library that is used to map the rates given in units of phisical time, to rates in units of computational Time-steps.
 * @param tipo the number that identifies the specie. 
 * */
void EscalaTiempoMetabolico(int tipo);

/** allocate the necesary memory for a lattice (system). 
 * @param[out] *es pointer to the lattice (system)
 * @param NDX the size of one side of the lattice: i=X coordinate 
 * @param NDY the size of the other side of the lattice: j=Y coordinate
 * @see estado 
 * */
void AlojaMemoria(estado *es, int NDX, int NDY);

/** Initialize a lattice (system) to be empty. 
 * @param *es pointer to the lattice (system).
 * */
void ResetEstado(estado *es);

/** Writtes to a <estado> a uniform random state of occupied sites of a species. 
 * @param *es lattice (system) where the random state will be "written".
 * @param[in] frac the fraction of the lattice that will be occupied for the generated random state.
 * @param[in] tipo the specie of the generated random state.
 * */
void GeneraEstadoAleatorio(estado *es, float frac, int tipo);

/**
 * Update an individual of a system <*es>. 
 * @param *es the system (lattice)
 * @param N the index of the individual that is going to be updated @see estado::INDICE
 * @param campo
 */
void ActualizaRyC(estado *es, int N, int campo);

/**
 * Makes the system <es> evolve one Time-step. 
 * One Time-step is counted when on average all the individuals in the system has been updated once.
 * @param *es pointer to the struct 'estado' that is going to be evolved. 
 */
void BarrMCcRyCamp(estado *es);
void EligeUniforme(int i,int j,int radio, sitio *vecino);
void InsertaIndividuosAleatorio(estado *es, int N, int tipo);
float OnPromRadio(estado *es, int radio);
float FuncionCorrelacion(estado *es,int radio);			//No eficiente
float FuncionCorrelacion2(estado *es,int radio);		//No eficiente
float Correlacion(estado *es,int radio);
float FuncionCorrelacionSpecies(estado *es,int radio,int TipoOrigen, int TipoDistante); //No eficiente
float CorrelacionEspecies(estado *es,int radio,int TipoOrigen, int TipoDistante);
int CuentaEspecie(estado *es, int tipo);
void InsertaIndividuoEn(estado *es,int i,int j,int tipo);
void AlojaMemoriaEspecie(int tipo);

void ActualizaRhoVsT_MP(estado *e,Float2D_MP *RhoVsT,Dist_MP *RhoDist);	
void IniciaMemoriaFloat2D_MP(Float2D_MP *ARRAY);
void IniciaMemoriaInt2D_MP(Int2D_MP *ARRAY);
void IniciaMemoriaDist_MP(Dist_MP *Dist);
void SumaFloat2D_MP(Float2D_MP *Origen, Float2D_MP *Destino);
void SumaDist_MP(Dist_MP *Origen, Dist_MP *Destino);
void InicializaFloat2D_MP(Float2D_MP *Objeto, int i_max, int j_max, int NoEnsambles);
void InicializaDist_MP(Dist_MP *Objeto, float TamParticion, float xIni, float xFin);
void SetSpecie(int NoEspecie, float Birth, float Coagulation, float Dead, float RadioBirth, float RadioCoa);
void ResetFloat2D_MP(Float2D_MP *ARRAY);
void ResetDist_MP(Dist_MP *Dist);
void SetSpecie2(int NoEspecie, float Birth, float Coagulation, float CoagulationIntra, float Dead, float RadioBirth, float RadioCoa, float RadioCoaIntra);
void ActualizaCorrelacion_MP(estado *es, Float1D_MP *corr);
void ActualizaCorrelacionTipo_MP(estado *es, Float1D_MP *corr, int TipoOrigen, int TipoObjetivo);
void SumaFloat1D_MP(Float1D_MP *Origen,Float1D_MP *Destino);
void InicializaFloat1D_MP(Float1D_MP *Objeto, int i_max);

void CorrXY(estado *es);
void CFFT(estado *es, Float2D_MP *correlacion);
void CFFT_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion);
void CompactaCorrelacion(Float2D_MP *corr2D, Float1D_MP *corrRadial);
void CFFT_Tipos_MP(estado *es, int NoEnsambles, Float2D_MP *correlacion,int TipoOrigen, int TipoDestino);
void ResetFloat1D_MP(Float1D_MP *ARRAY);
void DoblaCorrelacion(Float2D_MP *corr2D);
float Integra(Float1D_MP *Funcion, int inicial, int final);

void LiberaMemoria(estado *es);
void LiberaMemoriaFloat2D_MP(Float2D_MP *ARRAY);
void LiberaMemoriaFloat1D_MP(Float1D_MP *Objeto);
