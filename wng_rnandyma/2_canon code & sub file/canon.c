#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

// Row (OR) Column Communication.
void circular_shift(int *X, int* X_t, int N, MPI_Comm* comm, int size, int step, int pid)
{
	MPI_Request reqs[2];
	MPI_Status stats[2];

	int prev, next;
	int tag = 1;

	prev = (size + pid - step)%size;
	next = (size + pid + step)%size;
	//Communication between processors for the left shift and upward shift of Mat A and B respectively
	MPI_Irecv(X_t, (N*N), MPI_INT, next, tag, (*comm), &reqs[0]);
	MPI_Isend(X,   (N*N), MPI_INT, prev, tag, (*comm), &reqs[1]);

	MPI_Waitall(2, reqs, stats);

	return;
}



// Matrix Multiply and Accumulate.
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
			}
		}
	}
}

//Initialize Matrix.
void initMatrix(int* A, int N)
{
	int i, j;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			
			A[(N*i+j)] = 1;
			
		}
	}

	return;
}

// Print Matrix Elements.
void printMatrix(int* A, int N, char* str, int pid)
{
	int i, j;

	printf("%s (%02d)\n", str, pid);
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			printf("%d ", A[(N*i+j)]);
			
		}
		printf("\n");
	}
	printf("\n\n");
	return;
}



int main(int argc, char* argv[])
{
   int rpid, row_size, cpid, column_size, N, M;
	int pid, p;
   int i, j, k;

   int *A, *B, *C;
   int *Atemp, *Btemp;
   double time1;

   int Row_key, Column_key;
   MPI_Comm Row_comm, Column_comm;

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &pid);
   MPI_Comm_size(MPI_COMM_WORLD, &p);


   int iter = (int)sqrt(p); // every process should have m/sqrt(p) dim matrix
	

   if(argc == 2)
   {
      if(pid == 0) printf("\n\nnumP=%d, N=%s\n\n", p, argv[1]);
      M = atoi(argv[1]);
   }
   else
   {
      if(pid == 0)
      {
         printf("Command Line : %s", argv[0]);
         for(i = 1; i < argc; i++) printf(" %s", argv[i]);
      }
      MPI_Finalize();
		return 0;
   }
	N = M/iter;  // setting the size of the sub blocks
	A = (int*) malloc(N*N*sizeof(int)); // allocating  sub matrix block A for each proc
   	B = (int*) malloc(N*N*sizeof(int));// allocating sub matrix block B for each proc
   	C = (int*) malloc(N*N*sizeof(int)); // allocating matrix sub block C for each proc where C=AB
	Atemp = (int*) malloc(N*N*sizeof(int));// allocating Atemp matrix of same size of A to store the left shifted matrix
   	Btemp = (int*) malloc(N*N*sizeof(int)); // allocating Btemp matrix of same size of B to store the upward shifted matrix

	initMatrix(A, N);
	initMatrix(B, N);

	MPI_Datatype block, blocktype;
	int disp[p]; //required for gathering the final matrix subblocks
	int rcount[p]; //required for gathering the final matrix subblocks

	disp[0] = 0;
	rcount[0] = 1;
	for(i = 1; i < p; i++)
	{
		if((i%iter) == 0)
		{
			disp[i] = disp[(i-iter)] + (iter*N);
		}
		else
		{
			disp[i] = disp[(i-1)] + 1;
		}
		rcount[i] = 1;
	}
	MPI_Type_vector(N, N, M, MPI_INT, &block);
	MPI_Type_commit(&block);
	MPI_Type_create_resized(block, 0, (N*sizeof(int)), &blocktype);
	MPI_Type_commit(&blocktype);

	int *Z = NULL; // Final matrix to store all the computed values, only created in pid == 0
	if(pid == 0)
	{
       	    Z = (int*) malloc(M*M*sizeof(int));
	}


	time1 = MPI_Wtime();

	// Design Row Communicators.
	//R_color = pid/iter;
	Row_key = pid/iter;
	MPI_Comm_split(MPI_COMM_WORLD, Row_key, pid, &Row_comm);
	MPI_Comm_rank(Row_comm, &rpid);
	MPI_Comm_size(Row_comm, &row_size);

	// Design Column Communicators.
	//C_color = pid%iter;
	Column_key = pid%iter;
	MPI_Comm_split(MPI_COMM_WORLD, Column_key, pid, &Column_comm);
	MPI_Comm_rank(Column_comm, &cpid);
	MPI_Comm_size(Column_comm, &column_size);


	// Row Align A.
	circular_shift(A, Atemp, N, &Row_comm, row_size, Row_key, rpid);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Col Align B.
	circular_shift(B, Btemp, N, &Column_comm, column_size, Column_key, cpid);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Zero Initialize Once the Result Matrix.
	for(i = 0; i < (N*N); i++) C[i] = 0;

	for(i = 0; i < (N*N); i++)
	{
		A[i] = Atemp[i];
		B[i] = Btemp[i];
	}
	MatrixMul(A, B, C, N); // Continous Partial Product accumulator.
	
	for(k = 1; k < iter; k++)
	{
	   // Row Circular Shift Left A.
	   circular_shift(A, Atemp, N, &Row_comm, row_size, 1, rpid);
	   
	   
		// Col Circular Shift Up B.
	   circular_shift(B, Btemp, N, &Column_comm, column_size, 1, cpid);
	   
	   
		// Continous Partial Product accumulator.
		for(i = 0; i < (N*N); i++)
	   {
	   	A[i] = Atemp[i];
	   	B[i] = Btemp[i];
	   }
	   MatrixMul(A, B, C, N);
	}
	
	MPI_Comm_free(&Row_comm);
	MPI_Comm_free(&Column_comm);
	
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = MPI_Wtime() - time1;
	
	MPI_Gatherv(C, (N*N), MPI_INT, Z, rcount, disp, blocktype, 0, MPI_COMM_WORLD); // Gather C subblocks.
		
	//if(pid == 0) printMatrix(Z, M, "C", pid);

	if(pid == 0) printf("%02d : Elapsed Time (PARALLEL) is %fs\n", pid, time1);
	

	MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0;
}
