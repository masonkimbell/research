#ifndef GSL_MATRIX_LIB_H
#define GSL_MATRIX_LIB_H
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
//PRINTING THE MATRIX
void print_gsl_matrix(gsl_matrix *a, int n1, int n2){
	char matrix[(n1-1)*100];
	sprintf(matrix, "===start of matrix===");
	int i, j;
	for(i = 0; i < n1; i++){
		for(j = 0; j < n2; j++){
			sprintf(matrix, "%g ", gsl_matrix_get(a, i, j));
		}
		sprintf(matrix, "\n");
	}
	sprintf(matrix, "===end of matrix===");
}

//PRINTING THE EIGENVALUES
void print_eigs(gsl_vector_complex *v, int n){
	int i;
	for(i = 0; i < n; i++){
		printf("%g + %gi\n", GSL_REAL(gsl_vector_complex_get(v, i)), GSL_IMAG(gsl_vector_complex_get(v, i)));
	}
}

//writing all of the output to fname
//right now this function just prints but I'll deal with that later
void print_output(gsl_matrix *a, gsl_vector_complex *v, int n, int matid, char* out){
	sprintf(out, "%s===start of matrix===\n", out);
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			sprintf(out, "%s %g ", out, gsl_matrix_get(a, i, j));
		}
		sprintf(out, "%s\n", out);
	}
	sprintf(out, "%s===end of matrix===\n\n", out);

	sprintf(out, "%s~~~eigenvalues for %d x %d matrix #%d~~~\n", out, n, n, matid);
	for(i = 0; i < n; i++){
		sprintf(out, "%s%g + %gi\n", out, GSL_REAL(gsl_vector_complex_get(v, i)), GSL_IMAG(gsl_vector_complex_get(v, i)));
	}
	sprintf(out, "%s~~~end~~~\n\n", out);

	//printf("%s\n ^ has strlen %d\n", out, strlen(out));
}

//CREATING THE MATRIX
void gsl_offdiag(gsl_matrix *a, int n){
  int i;
  for(i = 0; i < n-1; i++){
    gsl_matrix_set(a, i, i+1, 1);
  }
}

void gsl_branch(gsl_matrix *a, int bp, int n){
	//bp = the leftmost point of the branch
	//three cases: top, bottom, and middle branch
	if(bp == -1)
		return;

	if(bp == 0){//top branch
		puts("top branch case");
		gsl_matrix_set(a, 0, 1, 0);
		gsl_matrix_set(a, 0, 2, 1);
	}else if(bp == n-2){//bottom branch
		puts("bottom branch case");
		gsl_matrix_set(a, n-3, n-1, 1);
		gsl_matrix_set(a, n-2, n-1, 0);
	}else if(bp > 0 && bp < n-2){//middle branch
		puts("middle branch case");
		int up = bp+1;
		int down = bp-1;
		gsl_matrix_set(a, down, up, 1); //for bp = 1: affects 0, 2
		gsl_matrix_set(a, bp, up, 0); //for bp = 1: affects 1, 2
		gsl_matrix_set(a, bp, up+1, 1); //for bp = 1: affects 1, 3
	}
}

//returns 1 if a zero is not found or 0 if one is found
int gsl_nonzero(gsl_matrix *a, int n){
  int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if(gsl_matrix_get(a, i, j) == 0){
        return 0;
      }
    }
  }
  return 1;
}

//MATRIX OPS

//add two square matrices
void gsl_add(gsl_matrix *A, gsl_matrix *B, gsl_matrix *result, int n){
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			gsl_matrix_set(result, i, j, (gsl_matrix_get(A, i, j) + gsl_matrix_get(B, i, j)));
		}
	}
}

//multiply two square matrices
void gsl_mult(gsl_matrix *A, gsl_matrix *B, gsl_matrix *result, int n){
	int i, j, k;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < n; k++){
				gsl_matrix_set(result, i, j, (gsl_matrix_get(result, i, j) + gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j)));
			}
		}
	}
}

//find the eigenvalues of a square matrix
void gsl_eigs(gsl_matrix *a, int size, int matrix_id, char* output){
	//hold the original matrix to write to file later
	gsl_matrix *a_copy;
	a_copy = gsl_matrix_calloc(size, size);
	gsl_matrix_memcpy(a_copy, a);

	//create workspace for computing eigenvalues
	gsl_eigen_nonsymm_workspace* w;
	w = gsl_eigen_nonsymm_alloc(size);
	gsl_eigen_nonsymm_params(0,0,w);

	//print_gsl_matrix(a, size, size);

	//create vector to hold the eigenvalues
	gsl_vector_complex *result;
	result = gsl_vector_complex_calloc(size);

	//calculate eigenvalues
	if(gsl_eigen_nonsymm(a, result, w) != 0){
		puts("error");
		exit(1);
	}

	//print eigenvalues
	//TODO: do something different with output
	//either return it or print it to a file
	
	print_output(a_copy, result, size, matrix_id, output);

	gsl_matrix_free(a_copy);
}

//SUBROUTINES FOR GCD
int max(int* arr, int n){
	int i;
	int max = arr[0];
	for(i = 1; i < n; i++){
		if(max < arr[i])
			max = arr[i];
	}
	return max;
}

int allZero(int* arr, int n){
	int i;
	for(i = 0; i < n; i++){
		if(arr[i] != 0)
			return 0;
	}
	return 1;
}

//find the gcd of the cycle lengths of matrix A, using up arrow n
int gcd(gsl_matrix *A, int n, int bp){
	//calculate cycle length for each cycle in graph
	//to deal with branches, bp is the leftmost point of the branch, if it exists
	//if end point of up arrow is...
	//above branch, subtract 1 from cycle length
	//on left side of branch, subtract 1 from cycle length
	//on right side of branch, do nothing
	//if start point of up arrow is...
	//on left side of branch, do nothing
	//on right side of branch, subtract 1 from cycle length
	//print_gsl_matrix(A, n, n);
	//printf("calculating cycle lengths of graph size %d...\n", n);
	int h, i, j;
	int len = 0;
	int* cycleLengths;
	cycleLengths = (int *) malloc(1 * sizeof(int));
	int count = 0;
		for(i = 1; i < n; i++){
			for(j = 0; j < i; j++){
			//	printf("i = %d and j = %d\n", i, j);
				if(gsl_matrix_get(A, i, j) == 1){
					len = i - j + 1;

					//do all checks for branches here
					//i = start point of up arrow
					//j = end point of up arrow
					if(bp > 0){
						if(j < bp && i > bp + 1)
							len--;
						else if(j == bp && i > bp + 1)
							len--;
						else if(j == bp + 1 && i > bp + 1)
							;//do nothing
	
						if(i == bp && j < bp)
							;//do nothing
						else if(i == bp + 1 && j < bp)
							len--;
					}

					count++;
					cycleLengths = (int *) realloc(cycleLengths, count * sizeof(int));
					cycleLengths[count - 1] = len;
				}
			}
		}
		//printf("found %d cycles with lengths...\n", count);

		int q;
		for(q = 0; q < count; q++){
			//printf("%d ", cycleLengths[q]);
		}

		//actual gcd calculation
		int gcd;
		int val = max(cycleLengths, count);
		while(val > 0){
			int* copy;
			copy = (int *) malloc(count * sizeof(int));
			int x;
			for(x = 0; x < count; x++){
				copy[x] = cycleLengths[x] % val;	
			}
			if(allZero(copy, count)){
				gcd = val;
				break;
			}
			free(copy);
			val--;
		}
		
		if(count == 0)
			gcd = 0;
		return gcd;
		free(cycleLengths);
}
#endif
