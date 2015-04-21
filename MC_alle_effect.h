#include "libPP_5.1.h"

typedef struct ALLE_ENV_STRUCT alle_env;

struct ALLE_ENV_STRUCT {
especie *param;
float Max_Metabolic;
int NumberOfSpecies;
float alle_effect;
};

void MC_sweep_alle(estado *es, alle_env *env);
void Update_alle(estado *es, alle_env *env, int N);
float calculate_metabolic_time(alle_env *env);
int store_sim_db(char *values);
