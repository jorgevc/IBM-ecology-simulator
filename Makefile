alle_effect: AlleEffect.c libPP_5.0.c EntSalArb_MP.c GNA.c conn_mysql.c
	gcc -o3 -fopenmp AlleEffect.c libPP_5.0.c EntSalArb_MP.c GNA.c conn_mysql.c -lfftw3 -lm `mysql_config --cflags --libs` -o alle.out
