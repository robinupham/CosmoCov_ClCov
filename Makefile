opts := -g -lgsl -lgslcblas -std=c99 -lm -ICosmoCov -funroll-loops -ffast-math -O3 -lfftw3

shear_clcov:
	gcc -o get_shear_clcov get_shear_clcov.c $(opts)
