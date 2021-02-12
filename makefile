# Makefile for multiple source files
# a.out: compiles research.c
# gsl_test: compiles gsl_test.c

# Make executable file called proj3
a.out: research.c
	mpicc research.c -lm

# Compile research code with GSL library
run_gsl: main.c
	mpicc -o run_gsl main.c -lgsl -lgslcblas -lm 
