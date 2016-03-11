Competencia: Competencia.c libPP_6.1.c EntSalArb_MP_Comp.c GNA.c
	gcc -fopenmp -O3 Competencia.c libPP_6.1.c EntSalArb_MP_Comp.c GNA.c -lfftw3 -lm -o crece.out

MCMC: MCMC.c MLE.c libPP_6.1.c EntSalArb_MP_Comp.c GNA.c model.c
	gcc -fopenmp -O3 MCMC.c libPP_6.1.c EntSalArb_MP_Comp.c GNA.c -lfftw3 -lm -o mcmc.out

