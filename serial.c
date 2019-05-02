
// C program to multiply two square matrices. 
#include <stdio.h> 
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
  
// This function multiplies mat1[][] and mat2[][], 
// and stores the result in res[][] 
void MatrixMul(int* A, int* B, int* C, int N)
{
	int i, j, k;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			for(k = 0; k < N; k++)
			{
				C[(N*i+j)] += A[(N*i+k)]*B[(N*k+j)];
			//	C[i][j] = A[i][k]*B[k][j];
			}
		}
	}
}
void initMatrix(int* A, int N)
{
	int i, j;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			//A[(N*i+j)] = rand()%3;
			A[(N*i+j)] = 1;
			//A[i][j] = 1;
		}
	}

	return;
}
  
int main(int argc, char* argv[]) 
{ 
    int i, j, k;
    int M;
    int pid,p;

   int *A, *B, *C;
   double time1, time2;
   //clock_t time1;
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &pid);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   //
   if(argc == 2)
   {
      printf("N=%s\n\n",argv[1]);
      M = atoi(argv[1]);
   }
	A = (int*) malloc(M*M*sizeof(int)); 
   	B = (int*) malloc(M*M*sizeof(int));
   	C = (int*) malloc(M*M*sizeof(int)); 
    initMatrix(A,M);
    initMatrix(B,M);
    //time1 = clock();
    time1 = MPI_Wtime();
    MatrixMul(A,B,C,M);
    time1 = MPI_Wtime()-time1;
    //time1 = clock()-time1;
    //double time_taken = ((double)time1)/CLOCKS_PER_SEC;
    printf("matrix mul took %f seconds\n", time1);
    /*printf("Result matrix is \n"); 
    for (i = 0; i < N; i++) 
    { 
        for (j = 0; j < N; j++) 
           printf("%d ", C[i][j]); 
        printf("\n"); 
    } */
  MPI_Finalize();
    return 0; 
} 

