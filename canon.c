#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
//#include <mkl.h>

// Row (OR) Column Communication.
void commMPI(int *X, int* X_t, int N, MPI_Comm* comm, int size, int step, int pid)
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


// Matrix Multiply
/*void mulMatrix(int* A, int* B, int* C, int N)
{
	int i, j, k;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			C[(N*i+j)] = 0;
  			//C[i][j] = 0;
			for(k = 0; k < N; k++)
			{
				C[(N*i+j)] += A[(N*i+k)]*B[(N*k+j)];
				//C[i][j] = A[i][k]*B[k][j];
			}
		}
	}
}*/

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
			//	C[i][j] = A[i][k]*B[k][j];
			}
		}
	}
}

// Random Initialize Matrix.
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
			//printf("%d ", A[i][j]);
			
		}
		printf("\n");
	}
	printf("\n\n");
	return;
}



int main(int argc, char* argv[])
{
   int rpid, rsize, cpid, csize, N, M;
	int pid, p;
   int i, j, k;

   int *A, *B, *C;
   int *Atemp, *Btemp;
   double time1, time2;

   int R_color, C_color;
   MPI_Comm R_comm, C_comm;

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
	A = (int*) malloc(N*N*sizeof(int)); // allocating matrix A for each proc
   	B = (int*) malloc(N*N*sizeof(int));// allocating matrix B for each proc
   	C = (int*) malloc(N*N*sizeof(int)); // allocating matrix C for each proc where C=AB
	Atemp = (int*) malloc(N*N*sizeof(int));// allocating Atemp matrix of same size of A to store the left shifted matrix
   	Btemp = (int*) malloc(N*N*sizeof(int)); // allocating Btemp matrix of same size of B to store the upward shifted matrix

	initMatrix(A, N);
	initMatrix(B, N);

	MPI_Datatype block, blocktype;
	int disp[p];
	int rcount[p];

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

	//int *X = NULL;
	//int *Y = NULL;
	int *Z = NULL;
	if(pid == 0)
	{
	   // X = (int*) malloc(p*N*N*sizeof(int));
       	   // Y = (int*) malloc(p*N*N*sizeof(int));
       	    Z = (int*) malloc(M*M*sizeof(int));
	}

/*
	MPI_Gatherv(A, (N*N), MPI_INT, X, rcount, disp, blocktype, 0, MPI_COMM_WORLD); // Gather A subblocks to matrix X for serial mul
	MPI_Gatherv(B, (N*N), MPI_INT, Y, rcount, disp, blocktype, 0, MPI_COMM_WORLD); // Gather B subblocks to matrix Y for serial mul

	if(pid == 0)
	{
	   time2 = MPI_Wtime();
		
		mulMatrix(X, Y, Z, (iter*N)); // Compute Result Matrix Sequentially for Verification.
	
		time2 = MPI_Wtime() - time2;
		printf("%02d : Elapsed Time (SERIAL) is   %fs\n", pid, time2);
		
		printMatrix(X, (iter*N), "X", pid);
	   printMatrix(Y, (iter*N), "Y", pid);
	   printMatrix(Z, (iter*N), "Z", pid);
		free(X);
		free(Y);
	}


	MPI_Barrier(MPI_COMM_WORLD);*/
	time1 = MPI_Wtime();

	// Design Row Communicators.
	R_color = pid/iter;
	MPI_Comm_split(MPI_COMM_WORLD, R_color, pid, &R_comm);
	MPI_Comm_rank(R_comm, &rpid);
	MPI_Comm_size(R_comm, &rsize);

	// Design Column Communicators.
	C_color = pid%iter;
	MPI_Comm_split(MPI_COMM_WORLD, C_color, pid, &C_comm);
	MPI_Comm_rank(C_comm, &cpid);
	MPI_Comm_size(C_comm, &csize);


	// Row Align A.
	commMPI(A, Atemp, N, &R_comm, rsize, R_color, rpid);
	MPI_Barrier(MPI_COMM_WORLD);
	
	// Col Align B.
	commMPI(B, Btemp, N, &C_comm, csize, C_color, cpid);
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
	   commMPI(A, Atemp, N, &R_comm, rsize, 1, rpid);
	   
	   
		// Col Circular Shift Up B.
	   commMPI(B, Btemp, N, &C_comm, csize, 1, cpid);
	   
	   
		// Continous Partial Product accumulator.
		for(i = 0; i < (N*N); i++)
	   {
	   	A[i] = Atemp[i];
	   	B[i] = Btemp[i];
	   }
	   MatrixMul(A, B, C, N);
	}
	
	MPI_Comm_free(&R_comm);
	MPI_Comm_free(&C_comm);
	
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = MPI_Wtime() - time1;
	
	MPI_Gatherv(C, (N*N), MPI_INT, Z, rcount, disp, blocktype, 0, MPI_COMM_WORLD); // Gather C subblocks.
	//MPI_Gather(C, N*N, MPI_INT, Z, N*N,MPI_INT, 0, MPI_COMM_WORLD);	
	//if(pid == 0) printMatrix(Z, M, "C", pid);

	if(pid == 0) printf("%02d : Elapsed Time (PARALLEL) is %fs\n", pid, time1);
	//if(pid == 0) printf("%02d : SPEEDUP (Serial vs Parallel) is %f\n", pid, time2/time1);

	MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0;
}
