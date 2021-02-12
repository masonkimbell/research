#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "matrixLib.c"
struct permutation{
  int* arr;
};

<<<<<<< HEAD
/*
void spewAdj(
*/
=======
>>>>>>> 8ed592e2d1c924a68ce24d58598570825748fe8c
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

//returns 1 if a zero is not found or 0 if one is found
int checkNonZero(floatMatrix m, int n){
  int i, j;
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      if(ACCESS2(m, i, j) == 0){
        return 0;
      }
    }
  }
  return 1;
}

int combine_arrays(int* arr1, int arr1sz, int* arr2, int arr2sz, int* result){  
  int i, j;
  for(i = 0; i < arr1sz; i++){
    result[i] = arr1[i];
  }

  for(j = 0; j < arr2sz; j++){
    result[j+arr1sz] = arr2[j];
  }
}

int max(int* arr, int n){
<<<<<<< HEAD

	int i;
	int max = arr[0];
	for(i = 1; i < n; i++){
		if(max < arr[i]) max = arr[i];
	}
	return max;
=======
  int i;
  int max = arr[0];
  for(i = 1; i < n; i++){
    if(max < arr[i]) max = arr[i];
  }
  return max;
>>>>>>> 8ed592e2d1c924a68ce24d58598570825748fe8c
}

int allZero(int* arr, int n){
  int i;
  for(i = 0; i < n; i++){
    if(arr[i] != 0)
      return 0;
  }
  return 1;
}

void gcd(floatMatrix A, int n, int bp){
<<<<<<< HEAD
	 


	// MPI File addition. 
	
	printf("HEREEEEE!!!!!\n");
	//calculate cycle length for each cycle in graph
	//to deal with branches, bp is the leftmost point of the branch, if it exists
	//if end point of up arrow is...
	//above branch, subtract 1 from cycle length
	//on left side of branch, subtract 1 from cycle length
	//on right side of branch, do nothing
	//if start point of up arrow is...
	//on left side of branch, do nothing
	//on right side of branch, subtract 1 from cycle length
	printFloatMatrix(&A);
	printf("calculating cycle lengths of graph size %d...\n", n);
	int h, i, j;
	int len = 0;
	int* cycleLengths;
	cycleLengths = (int *) malloc(1 * sizeof(int));
	int count = 0;
	//for(h = 0; h < n-1; h++){
	for(i = 1; i < n; i++){
		for(j = 0; j < i; j++){
			//	printf("i = %d and j = %d\n", i, j);
			if(ACCESS2(A, i, j) == 1){
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
	printf("found %d cycles with lengths...\n", count);
	print_array(cycleLengths, count);

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
	printf("the gcd is %d\n", gcd);

	free(cycleLengths);
=======
  //calculate cycle length for each cycle in graph
  //to deal with branches, bp is the leftmost point of the branch, if it exists
  //if end point of up arrow is...
  //above branch, subtract 1 from cycle length
  //on left side of branch, subtract 1 from cycle length
  //on right side of branch, do nothing
  //if start point of up arrow is...
  //on left side of branch, do nothing
  //on right side of branch, subtract 1 from cycle length
  printFloatMatrix(&A);
  printf("calculating cycle lengths of graph size %d...\n", n);
  int h, i, j;
  int len = 0;
  int* cycleLengths;
  cycleLengths = (int *) malloc(1 * sizeof(int));
  int count = 0;
  //for(h = 0; h < n-1; h++){
  for(i = 1; i < n; i++){
    for(j = 0; j < i; j++){
      //	printf("i = %d and j = %d\n", i, j);
      if(ACCESS2(A, i, j) == 1){
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
  printf("found %d cycles with lengths...\n", count);
  print_array(cycleLengths, count);

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
  printf("the gcd is %d\n", gcd);

  free(cycleLengths);
>>>>>>> 8ed592e2d1c924a68ce24d58598570825748fe8c
}

//run instructions: ./a.out n
int main(int argc, char** argv){
<<<<<<< HEAD
	MPI_Init(NULL, NULL);

	if( argc < 2 ){
		puts("Usage: ./a.out N");
		return 1;
	}

	// TODO: make sure world size is power of 2!!
  // File Output of adj. Matricies. 

  int err = 0; 
	
	
	int world_size;

  // *** File Handling MPI. 
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Comm_size(world, &world_size);
	MPI_Info infoin, infout;
	MPI_File fh;
  int me;
  MPI_Comm_rank(world, &me);


	// Error with MPI File Open function. 

 MPI_File_open(world, "adjTest.txt", MPI_MODE_RDWR | MPI_MODE_CREATE, infoin, &fh);
 //if(err)  MPI_Abort( MPI_COMM_WORLD, 911 );



	printf("Starting node %d in world size: %d\n", me, world_size);

	//n = size of matrix
	int n = atoi(argv[1]);
	//loop var
	int i;

	int nonzero = 0;
	int total = 0;

	//k = number of entries in lower triangular portion
	int k = n*(n-1)/2;

	//world size is 8, so the first 3 values are fixed since 2^3 = 8
	int prefix_size = (int) (log2(world_size));
	k -= prefix_size;

	int count = 0;

	//arr: holds a binary permutation for lower tri
	int* suffix; // suffix array
	suffix = (int*) malloc(k*sizeof(int));

	//matrix: holds graph info
	floatMatrix A;
	//create an nxn matrix of all 0s
	initFloatMatrix(&A, n, n, 1);//all zeros

	// set top right to 1's
	makeOffDiagonal(&A, n);

	int branchPoint = 2;
	//^ maybe make this a parameter?

	/* FOR TESTING THAT createBranch WORKS PROPERLY
	*/createBranch(&A, branchPoint, n);
	puts("BRANCH");
	printFloatMatrix(&A);
	puts("END BRANCH");
	/*END TESTING */

	//TODO: use a new parameter to pick what branch we look at
	//change the above line based on the value of that parameter

	int possibilities = (1 << k);

	//printf("me = %d\n", me);
	int prefix[prefix_size];
	convert_to_bin(prefix, me, prefix_size);
	//print_array(firstThree, 3);

	int w = 0;

	//
	//
	// parallel stuff starts here
	//
	//

	//creates matrices with each lower tri possibility
	for(i = 0; i < possibilities; i++){
		convert_to_bin(suffix, i, k);

		int *combined;
		combined = (int*) malloc((k+prefix_size)*sizeof(int));
		combine_arrays(prefix, prefix_size, suffix, k, combined);

		if(i%100 == 0){
			printf("Node %d Starting possibility:",me);
			print_array(combined,k+prefix_size);
			puts("");
		}

		int a, b;
		int count = 0;
		for(a = 1; a < n; a++){
			for(b = 0; b < a; b++){
				//printf("adding %d, %d\n", a, b);
				ACCESS2(A, a, b) = combined[count]; 
				count++;
			}
		}

		//prod: the product matrix
		floatMatrix prod;
		initFloatMatrix(&prod, n, n, 1);
		copyFloatMatrix(&prod, &A);

		//sum: holds A^1 + A^2 + ... + A^n
		//initially starts as A^1
		floatMatrix sum;
		initFloatMatrix(&sum, n, n, 1);
		copyFloatMatrix(&sum, &A);
		int pow = 1;
		int x,y;
		int numMult = 1;

		//calculate A^n
		//need to accumulate after each time (A^1 + A^2 + ... + A^n)
		//outer loop: sum
		//inner loop: product

		//confirm that A holds the right matrix (it does)

		if (me == 3){
			printf("\nSTARTING\n");
			printFloatMatrix(&A);
			printf("MATRIX");
		}

		/* FOR TESTING THAT createBranch WORKS PROPERLY
			 createBranch(&A, 2, n);
			 puts("BRANCH");
			 printFloatMatrix(&A);
			 puts("END BRANCH");
			 END TESTING */


		//counters to keep track of # of strongly connected
		gcd(A, n, branchPoint);

		for(x = 1; x <= n; x++){
			floatMatrix result;
			initFloatMatrix(&result, n, n, 1);
			int numberofmult = 1;


			//compute A^y
			for(y = 1; y <= x; y++){
				//we need to keep A, accumulate product in prod
				//prod = A * prod, where prod is the previous product
				if(y == 1){

					copyFloatMatrix(&prod, &A);

					}

				else{
					numberofmult++;
					mult_seq(&A, &prod, &result, n);
					copyFloatMatrix(&prod, &result);
				}
			}

			if(x!=1) add_seq(&sum, &prod, &sum, n);

			free(result.arr);
		}



		if(checkNonZero(sum, n)){
			nonzero++;
			puts("this graph is strongly connected");
		}
		//printf("out\n");
		total++;
	
		//free all of the matrices used: A, prod, result, sum
		free(prod.arr);
		free(sum.arr);
		free(combined);
	}


	gcd(A, n, branchPoint);
	printf("%d/%d are strongly connected\n", nonzero, total);
	printf("Here\n");
	free(A.arr);
	free(suffix);
	//todo:
	//put into parallel - in progress
	//test eigenvalue function - just add back eigenvalue stuff once parallelization works
	MPI_Finalize();

	return 0;
=======
  MPI_Init(NULL, NULL);

  if( argc < 2 ){
    puts("Usage: ./a.out N");
    return 1;
  }

  // TODO: make sure world size is power of 2!!

  int world_size;
  MPI_Comm world = MPI_COMM_WORLD;
  MPI_Comm_size(world, &world_size);

  int me;
  MPI_Comm_rank(world, &me);

  printf("Starting node %d in world size: %d\n", me, world_size);

  //n = size of matrix
  int n = atoi(argv[1]);
  //loop var
  int i;

  int nonzero = 0;
  int total = 0;

  //k = number of entries in lower triangular portion
  int k = n*(n-1)/2;

  //world size is 8, so the first 3 values are fixed since 2^3 = 8
  int prefix_size = (int) (log2(world_size));
  k -= prefix_size;

  int count = 0;

  //arr: holds a binary permutation for lower tri
  int* suffix; // suffix array
  suffix = (int*) malloc(k*sizeof(int));

  //matrix: holds graph info
  floatMatrix A;
  //create an nxn matrix of all 0s
  initFloatMatrix(&A, n, n, 1);//all zeros

  // set top right to 1's
  makeOffDiagonal(&A, n);

  int branchPoint = 2;
  //^ maybe make this a parameter?

  /* FOR TESTING THAT createBranch WORKS PROPERLY
  */createBranch(&A, branchPoint, n);
  puts("BRANCH");
  printFloatMatrix(&A);
  puts("END BRANCH");
  /*END TESTING */

  //TODO: use a new parameter to pick what branch we look at
  //change the above line based on the value of that parameter

  int possibilities = (1 << k);

  //printf("me = %d\n", me);
  int prefix[prefix_size];
  convert_to_bin(prefix, me, prefix_size);
  //print_array(firstThree, 3);

  int w = 0;

  //
  //
  // parallel stuff starts here
  //
  //

  //creates matrices with each lower tri possibility
  for(i = 0; i < possibilities; i++){
    convert_to_bin(suffix, i, k);

    int *combined;
    combined = (int*) malloc((k+prefix_size)*sizeof(int));
    combine_arrays(prefix, prefix_size, suffix, k, combined);

    if(i%100 == 0){
      printf("Node %d Starting possibility:",me);
      print_array(combined,k+prefix_size);
      puts("");
    }

    int a, b;
    int count = 0;
    for(a = 1; a < n; a++){
      for(b = 0; b < a; b++){
        //printf("adding %d, %d\n", a, b);
        ACCESS2(A, a, b) = combined[count]; 
        count++;
      }
    }
    /*
       if(me == 4){
       printf("================TESTING %d=====================\n", me);
       printFloatMatrix(&A);
       printf("================END=====================\n");
       }
       */
    //
    //
    // works up to this point
    //
    // think the issue has something to do with the "outer" world variable being used in the multiply function
    // when executed, no computation happens at all
    //
    //

    //prod: the product matrix
    floatMatrix prod;
    initFloatMatrix(&prod, n, n, 1);
    copyFloatMatrix(&prod, &A);

    //sum: holds A^1 + A^2 + ... + A^n
    //initially starts as A^1
    floatMatrix sum;
    initFloatMatrix(&sum, n, n, 1);
    copyFloatMatrix(&sum, &A);
    int pow = 1;
    int x,y;
    int numMult = 1;

    //calculate A^n
    //need to accumulate after each time (A^1 + A^2 + ... + A^n)
    //outer loop: sum
    //inner loop: product

    //confirm that A holds the right matrix (it does)

    if (me == 3){
      printf("\nSTARTING\n");
      printFloatMatrix(&A);
      printf("MATRIX");
    }

    /* FOR TESTING THAT createBranch WORKS PROPERLY
       createBranch(&A, 2, n);
       puts("BRANCH");
       printFloatMatrix(&A);
       puts("END BRANCH");
       END TESTING */


    //counters to keep track of # of strongly connected

    for(x = 1; x <= n; x++){
      floatMatrix result;
      initFloatMatrix(&result, n, n, 1);
      int numberofmult = 1;


      //compute A^y
      for(y = 1; y <= x; y++){
        //we need to keep A, accumulate product in prod
        //prod = A * prod, where prod is the previous product
        if(y == 1){
          copyFloatMatrix(&prod, &A);
          //printf("=====A^%d======\n",numberofmult);
          //printFloatMatrix(&prod);
          //printf("case 1\n");
        }else{
          numberofmult++;
          mult_seq(&A, &prod, &result, n);
          copyFloatMatrix(&prod, &result);
          //prod will hold A^x
          //printf("=====A^%d======\n",numberofmult);
          //printFloatMatrix(&prod);
          //printf("case 2\n");
        }
      }

      //printf("adding A^%d", numberofmult);
      //sum = prod + sum, where prod is the previous product
      if(x!=1) add_seq(&sum, &prod, &sum, n);

      //printf("=====sum=====\n");
      //printFloatMatrix(&sum);
      free(result.arr);
    }
    gcd(A, n, branchPoint);
    /*
       printf("=====product=====\n");
       printFloatMatrix(&prod);

       if(me == 4){
       printf("=====sum=====\n");
       printFloatMatrix(&sum);
       }
       */

    if(checkNonZero(sum, n)){
      nonzero++;
      puts("this graph is strongly connected");
    }
    //printf("out\n");
    total++;

    //free all of the matrices used: A, prod, result, sum
    free(prod.arr);
    free(sum.arr);
    free(combined);
  }

  printf("%d/%d are strongly connected\n", nonzero, total);
  free(A.arr);
  free(suffix);
  //todo:
  //put into parallel - in progress
  //test eigenvalue function - just add back eigenvalue stuff once parallelization works
  MPI_Finalize();

  return 0;
>>>>>>> 8ed592e2d1c924a68ce24d58598570825748fe8c
}
