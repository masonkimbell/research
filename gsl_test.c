#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

void print_gsl_matrix(gsl_matrix *a, int n1, int n2){
	puts("===start of matrix===");
	int i, j;
	for(i = 0; i < n1; i++){
		for(j = 0; j < n2; j++){
			printf("%g ", gsl_matrix_get(a, i, j));
		}
		puts("");
	}
	puts("===end of matrix===");
}

void print_gsl_vector_complex(gsl_vector_complex *v, int n){
	puts("===start of vector===");
	int i;
	for(i = 0; i < n; i++){
		printf("%g + %gi\n", GSL_REAL(gsl_vector_complex_get(v, i)), GSL_IMAG(gsl_vector_complex_get(v, i)));
	}
	puts("===end of vector===");
}

int main(){
	const int size = 4;

	//create workspace for computing eigenvalues
	gsl_eigen_nonsymm_workspace* w;
	w = gsl_eigen_nonsymm_alloc(size);
	gsl_eigen_nonsymm_params(0,0,w);

	//create the actual matrix
	gsl_matrix *a;
	a = gsl_matrix_calloc(size,size);
	
	//sets (0, 0) to 1, (1, 1) to 2, and (2, 2) to 1
	gsl_matrix_set(a, 0, 0, 4);
	gsl_matrix_set(a, 0, 1, 3);
	gsl_matrix_set(a, 0, 2, 2);
	gsl_matrix_set(a, 0, 3, 9);
	gsl_matrix_set(a, 1, 0, 8);
	gsl_matrix_set(a, 1, 1, 7);
	gsl_matrix_set(a, 1, 2, 2);
	gsl_matrix_set(a, 1, 3, 4);
	gsl_matrix_set(a, 2, 0, 1);
	gsl_matrix_set(a, 2, 1, 6);
	gsl_matrix_set(a, 2, 2, 3);
	gsl_matrix_set(a, 2, 3, 7);
	gsl_matrix_set(a, 3, 0, 6);
	gsl_matrix_set(a, 3, 1, 3);
	gsl_matrix_set(a, 3, 2, 1);
	gsl_matrix_set(a, 3, 3, 1);
	
	print_gsl_matrix(a, size, size);

	//create vector to hold the eigenvalues
	gsl_vector_complex *result;
	result = gsl_vector_complex_calloc(size);

	//calculate eigenvalues
	if(gsl_eigen_nonsymm(a, result, w) != 0){
		puts("error");
		return 1;
	}

	//print eigenvalues
	print_gsl_vector_complex(result, size);

	//note: maybe make a function to convert an old matrix to gsl?


	//free stuff
	gsl_matrix_free(a);
	return 0;
}
