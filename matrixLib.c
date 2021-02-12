#ifndef MATRIXLIB_C
#define MATRIXLIB_C

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <mpi.h>
#include <math.h>

#define INDEX(m,n,i,j) n*i + j
#define ACCESS(A,i,j) A->arr[INDEX(A->rows, A->cols, i, j)]
#define ACCESS2(A,i,j) A.arr[INDEX(A.rows, A.cols, i, j)]

typedef struct matrix{
	int rows, cols;
	int *arr;
}matrix;

typedef struct floatMatrix{
	int rows, cols;
	float *arr;
} floatMatrix;

typedef struct spot{
	int i, j;
}spot;

//print given by Dr. Anderson
void printMatrix(matrix *A){
	int i,j;
	for (i=0; i < A->rows; i++){
		for (j=0; j < A->cols; j++){
			printf("%d ", ACCESS(A,i,j));
		}
		puts("");
	}
}

void printFloatMatrix(floatMatrix *A){
	int i, j;
	for (i = 0; i < A->rows; i++)
	{
		for (j = 0; j < A->cols; j++)
		{
			printf("%1.0f ", ACCESS(A, i, j));
		}
		puts("");
	}
}

// sets all upper-triangle entries to 1
void makeOffDiagonal(floatMatrix *m, int n){
	int i;
	for(i = 0; i < n-1; i++){
		ACCESS(m, i, i+1) = 1;
	}
}

// create branch starting at point bp
void createBranch(floatMatrix *m, int bp, int n){
	//bp = the leftmost point of the branch
	//three cases: top, bottom, and middle branch
	if(bp == -1)
		return;

	if(bp == 0){//top branch
		puts("top branch case");
		ACCESS(m, 0, 1) = 0;
		ACCESS(m, 0, 2) = 1;
	}else if(bp == n-2){//bottom branch
		puts("bottom branch case");
		ACCESS(m, n-3, n-1) = 1;
		ACCESS(m, n-2, n-1) = 0;
	}else if(bp > 0 && bp < n-2){//middle branch
		puts("middle branch case");
		int up = bp+1;
		int down = bp-1;
		ACCESS(m, down, up) = 1; //for bp = 1: affects 0, 2
		ACCESS(m, bp, up) = 0; //for bp = 1: affects 1, 2
		ACCESS(m, bp, up+1) = 1; //for bp = 1: affects 1, 3
	}
}

void normalize(floatMatrix *p, int n)
{
	//float w = sqrt( p->x * p->x + p->y * p->y + p->z * p->z );

	int i, j, k;
	float w = 0;
	for (i=0; i < n; i++)
	{
		w += ACCESS(p, i, 0)*ACCESS(p, i, 0);
	}
	w = sqrt(w);
	//  printf("Value of w: %f\n", w);
	for (i=0; i < n; i++)
	{
		ACCESS(p, i, 0) /= w;
	}
	//printf("\n\n");
	//printFloatMatrix(p);
	//printf("\n");
}


//helper to convert a matrix entry to a i, j spot to then be used with ACCESS
spot convertToIJ(int place, int colSize)
{
	spot temp;
	temp.i = place/colSize;
	temp.j = place % colSize;
	return temp;
}
//init given by Dr. Anderson
void initMatrix(matrix *A, int r, int c, int empty)
{
	A->rows = r;
	A->cols = c;
	A->arr = malloc(r*c*sizeof(int));

	if (empty == 0)
	{
		int i, j;
		for (i=0; i < r; i++)
			for (j=0; j < c; j++)
				ACCESS(A,i,j) = rand() % 10 + 1;
	}
	else
	{
		int i, j;
		for (i=0; i < r; i++)
			for (j=0; j < c; j++)
				ACCESS(A,i,j) = 0;
	}
}
void initFloatMatrix(floatMatrix *A, int r, int c, int empty)
{
	A->rows = r;
	A->cols = c;
	A->arr = malloc(r * c * sizeof(float));

	if (empty == 0)
	{
		int i, j;
		for (i = 0; i < r; i++)
			for (j = 0; j < c; j++)
				ACCESS(A, i, j) =  (float) (rand() % 20 + 1);
	}
	else if (empty == 1)
	{
		int i, j;
		for (i = 0; i < r; i++)
			for (j = 0; j < c; j++)
				ACCESS(A, i, j) = 0.0;
	}
	else if (empty == 2)//identity
	{
		int i, j;
		for (i = 0; i < r; i++)
			for (j = 0; j < c; j++)
				if (i == j)
					ACCESS(A, i, j) = 1.0;
				else
					ACCESS(A, i, j) = 0.0;
	}
	else if (empty == 3)//for testing
	{
		ACCESS(A, 0, 0) = 0.0;
		ACCESS(A, 1, 0) = 1.0;
		ACCESS(A, 2, 0) = 0.0;
		ACCESS(A, 0, 1) = 0.0;
		ACCESS(A, 1, 1) = 0.0;
		ACCESS(A, 2, 1) = 0.0;
		ACCESS(A, 0, 2) = 0.0;
		ACCESS(A, 1, 2) = 2.0;
		ACCESS(A, 2, 2) = 2.0;
	}
	else if (empty == 4)
	{

		int i, j;
		for (i = 0; i < r; i++)
			for (j = 0; j < c; j++)
				ACCESS(A, i, j) = 1.0;
	}
}

void copyFloatMatrix(floatMatrix *B, floatMatrix *C)
{
	int i, j;
	for (i = 0; i < B->rows; i++)
		for (j = 0; j < B->cols; j++)
			ACCESS(B, i, j) =  ACCESS(C, i, j);
}

void add_seq(floatMatrix *A, floatMatrix *B, floatMatrix *C, int n){
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			ACCESS(C, i, j) = ACCESS(A, i, j) + ACCESS(B, i, j);
		}
	}
}

void mult_seq(floatMatrix *A, floatMatrix *B, floatMatrix *C, int n){
	int i, j, k;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			for(k = 0; k < n; k++){
				ACCESS(C, i, j) = ACCESS(C, i, j) + ACCESS(A, i, k) * ACCESS(B, k, j);
			}
		}
	}
}
#endif
