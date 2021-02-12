#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include "gsl_matrix_lib.h"

void print_array(int* arr, int n){
	int i;
  	for(i = 0; i < n; i++){
    		printf("%d ", arr[i]);
  	}
  	printf("\n");
}

void swap(int *x, int *y){
  	int temp;
  	temp = *x;
	*x = *y;
  	*y = temp;
}

//takes val and converts it to a binary array of size n called arr
void convert_to_bin(int* arr, int val, int n){
	int count;
  	int count2 = 0;
  	for(count = n-1; count >= 0; count--){
    		//printf("comparing %d, %d\n", val, (1<<count));
    		if((1 << count) <= val){
      			arr[count2] = 1;
      			val = val - (1 << count);
    		}else{
      			arr[count2] = 0;
    		}
    		count2++;
  	}
  	//printf("val = %d\n", val);
  	//print_array(arr, n);
}

void combine_arrays(int* arr1, int arr1sz, int* arr2, int arr2sz, int* result){  
  	int i, j;
  	for(i = 0; i < arr1sz; i++){
    		result[i] = arr1[i];
  	}

  	for(j = 0; j < arr2sz; j++){
    		result[j+arr1sz] = arr2[j];
  	}
}

//run instructions: ./a.out n
int main(int argc, char** argv){
  	MPI_Init(NULL, NULL);

  	if( argc < 2 ){
    		puts("Usage: ./a.out N");
    		exit(1);
 	}
	
  // TODO: make sure world size is power of 2!!

  	int world_size;
  	MPI_Comm world = MPI_COMM_WORLD;
  	MPI_Comm_size(world, &world_size);

  	int me;
  	MPI_Comm_rank(world, &me);
	
	MPI_Offset offset;

	offset = me * ((atoi(argv[1])-1) * 70 + 16);

	//make sure the offset is different for each node
	//printf("offset for node %d is %d\n", me, offset);

	MPI_File fh;
	MPI_Status status;

	char filename[20];
	sprintf(filename, "Out/outfile%d", atoi(argv[1]));
	if(MPI_File_open(world, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh) < 0){
		puts("error: need Out directory");
		exit(1);
	}
	//test that parallelization works
  	//printf("Starting node %d in world size: %d\n", me, world_size);

  	//n = size of matrix
  	long long n = atoi(argv[1]);
  	//loop var
	//will hold the number of matrix we are testing, the number is determined by the lower triangular portion
	//ex: lower triangular array 1 0 0 1 0 = 18
  	long long i;

  	long long nonzero = 0;
  	long long total = 0;

  	//k = number of entries in lower triangular portion
  	long long k = n*(n-1)/2;

  	//world size is 8, so the first 3 values are fixed since 2^3 = 8
  	int prefix_size = (int) (log2(world_size));
  	k -= prefix_size;

  	int count = 0;

  	//arr: holds a binary permutation for lower tri
  	int* suffix; // suffix array
  	suffix = (int*) malloc(k*sizeof(int));

	//everything above here works without the gsl stuff
	//just changing all matrices to gsl

  	//matrix: holds graph info
  	//create an nxn matrix of all 0s
  	gsl_matrix *A;
  	A = gsl_matrix_calloc(n,n);

  	//make matrix A off diagonal (set top right to 1's)
  	gsl_offdiag(A, n);

  	int branchPoint = 0;
  	//^ maybe make this a parameter?

  	//TODO: change this to work with gsl
  	/* FOR TESTING THAT createBranch WORKS PROPERLY
  	gsl_branch(A, branchPoint, n);
  	print_gsl_matrix(A, n, n);
  	END TESTING */

	long long one = 1;
  	long long possibilities = (long long)(one << k);
	printf("possibilities = %lld\n", possibilities);
  	//printf("me = %d\n", me);
  	int prefix[prefix_size];
  	convert_to_bin(prefix, me, prefix_size);
  	//print_array(firstThree, 3);

  	int w = 0;

  	//
  	// parallel stuff starts here
  	//

  	//creates matrices with each lower tri possibility
	for(i = 0; i < possibilities; i++){
		if(i == possibilities / 5 * 1)
			puts("20 percent done");
		else if(i == possibilities / 5 * 2)
			puts("40 percent done");
		else if(i == possibilities / 5 * 3)
			puts("60 percent done");
		else if(i == possibilities / 5 * 4)
			puts("80 percent done");
	
		//printf("TESTING MATRIX %d...\n", i);
    		convert_to_bin(suffix, i, k);
	    	int *combined;
    		combined = (int*) malloc((k+prefix_size)*sizeof(int));
    		combine_arrays(prefix, prefix_size, suffix, k, combined);

		/*
    		if(i%100 == 0){
      			printf("Node %d Starting possibility:",me);
      			print_array(combined,k+prefix_size);
      			puts("");
    		}
		*/

    		int a, b;
    		int count = 0;
    		for(a = 1; a < n; a++){
      			for(b = 0; b < a; b++){
        			//printf("adding %d, %d\n", a, b);
        			gsl_matrix_set(A, a, b, combined[count]);
        			count++;
      			}
    		}

    		//prod: the product matrix
    		gsl_matrix *prod;
    		prod = gsl_matrix_calloc(n, n);
   		gsl_matrix_memcpy(prod, A);

    		//sum: holds A^1 + A^2 + ... + A^n
    		//initially starts as A^1
    		gsl_matrix *sum;
    		sum = gsl_matrix_calloc(n, n);
    		gsl_matrix_memcpy(sum, A);
    
    		int pow = 1;
    		int x,y;
   	 	int numMult = 1;

    		//calculate A^n
    		//need to accumulate after each time (A^1 + A^2 + ... + A^n)
    		//outer loop: sum
    		//inner loop: product
		
		//make sure A holds the correct matrix
   		//print_gsl_matrix(A, n, n);    


    		//counters to keep track of # of strongly connected

    		for(x = 1; x <= n; x++){
	    		gsl_matrix *result;
	    		result = gsl_matrix_calloc(n, n);
  	    		int numberofmult = 1;


      			//compute A^y
    			for(y = 1; y <= x; y++){
        			//we need to keep A, accumulate product in prod
        			//prod = A * prod, where prod is the previous product
        			if(y == 1){
					gsl_matrix_memcpy(prod, A);
        	  			//print_gsl_matrix(prod, n, n);
       		  			//printf("case 1\n");
        			}else{
          				numberofmult++;
          				gsl_mult(A, prod, result, n);
          				gsl_matrix_memcpy(prod, result);
          				//prod will hold A^x

	     				//print_gsl_matrix(prod, n, n);
          				//printf("case 2\n");
        			}
      			}

     		 	//printf("adding A^%d", numberofmult);
      			//sum = prod + sum, where prod is the previous product
      			if(x!=1)
        			gsl_add(sum, prod, sum, n);
      			//printf("=====sum=====\n");
      			//print_gsl_matrix(sum, n, n);
      			gsl_matrix_free(result);
    		}
		
    		/*
        	printf("=====product=====\n");
        	print_gsl_matrix(prod, n, n);
        	*/

		//STRONGLY CONNECTED TEST
    		if(gsl_nonzero(sum, n)){
      			nonzero++;
      			//printf("Matrix %d is strongly connected\n", i);

			//only check eigenvalues of strongly connected matrices
			//TODO: check eigenvalues of strongly connected matrices where gcd = k
			//TODO: change output to file instead of print to terminal

			//find the gcd of strongly connected matrix
			int gcd_val;
			gcd_val = gcd(A, n, branchPoint);
			if(gcd_val > 2){
				//need to convert int to string and then combine them all

				int fsize = (n-1) * 70 + 16;
				char output[fsize];
				sprintf(output, "Proc = %d gcd = %d\n", me, gcd_val);

				//printf("SC Matrix w/ID#%d has gcd = %d\n", i, gcd_val);
				gsl_matrix *a_copy;
				a_copy = gsl_matrix_calloc(n,n);
				gsl_matrix_memcpy(a_copy, A);
				
				gsl_eigs(a_copy, n, i, output);
				
				MPI_File_write_at(fh, offset, output, fsize, MPI_CHAR, &status);
				offset += fsize * world_size;
				gsl_matrix_free(a_copy);
				//puts("");
			}
    		}
    		//printf("out\n");
    		total++;
    		//free all of the matrices used: A, prod, result, sum
    		gsl_matrix_free(prod);
    		gsl_matrix_free(sum);
    		free(combined);
  	}
	MPI_File_close(&fh);
	//test output for each individual processor
  	//printf("proc %d: %d/%d are strongly connected\n", me, nonzero, total);
	int total_nonzero;
	//combine output from all processors into one integer
	MPI_Reduce(&nonzero, &total_nonzero, 1, MPI_INT, MPI_SUM, 0, world);
	if(me == 0)
		printf("%d/%lld are strongly connected\n", total_nonzero, total*world_size);
  	gsl_matrix_free(A);
  	free(suffix);
  	MPI_Finalize();

}
