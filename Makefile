alle_effect: AlleEffect.c MC_alle_effect.c libPP_5.1.c EntSalArb_MP.c GNA.c
	gcc -o3 -fopenmp AlleEffect.c MC_alle_effect.c libPP_5.1.c EntSalArb_MP.c GNA.c -lfftw3 -lm -lsqlite3 -o alle.out
