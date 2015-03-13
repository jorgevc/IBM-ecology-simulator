alle_effect: AlleEffect.c MC_alle_effect.c libPP_5.1.c EntSalArb_MP.c GNA.c conn_mysql.c
	gcc -o3 -fopenmp AlleEffect.c MC_alle_effect.c libPP_5.1.c EntSalArb_MP.c GNA.c conn_mysql.c -lfftw3 -lm `mysql_config --cflags --libs` -o alle.out
