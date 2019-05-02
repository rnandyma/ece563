#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void commMPI(int *X, int* X_t, int N, MPI_Comm* comm, int p, int gpid)
{
	MPI_Request reqs[2];
	MPI_Status stats[2];

	int prev, next;
	int tag = 1;
        if(gpid == 0)
	prev = p -1;
	else
	prev = gpid-1;
        next = (gpid + 1)%p;
	//prev = (size + pid - step)%size;
	//next = (size + pid + step)%size;
	//Communication between processors for the left shift and upward shift of Mat A and B respectively
	MPI_Irecv(X_t, (N*N/p), MPI_INT, next, tag, (*comm), &reqs[0]);
	MPI_Isend(X,   (N*N/p), MPI_INT, prev, tag, (*comm), &reqs[1]);

	MPI_Waitall(2, reqs, stats);

	return;
}

// Multiply NxM matrixA with MxQ matrixB to get NxQ matrix C.
void mulMatrix(int *A, int *B, int *C, int N, int M, int shift, int width, int pid)
{
	int i, j, k;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			//C[(N*i+j)] = 0;
			C[(i*M)+(j+(shift*width)+(pid*width))%M] = 0;
			for(k = 0; k < M; k++)
			{
				//C[(N*i+j)] += A[(M*i+k)] * B[(N*k+j)];
				C[(i*M)+(j+(shift*width)+(pid*width))%M] += A[(M*i+k)] * B[(j*M+k)];
			}
		}
	}

	return;
}


void initSubMatrix(int *A, int N)
{
	int i;

	for(i = 0; i < N; i++)
	{
		//A[i] = rand()%3;
 		A[i] = i%4;
	}
	return;
}

void printMatrix(int *A, int N, int M, char *str, int gpid)
{
	int i, j;

	printf("%s (%02d)\n", str, gpid);
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < M; j++) printf("%d ", A[M*i+j]);
	   printf("\n");
	}

	return;
}


int main(int argc, char* argv[])
{
	int i, j, k;
	int gpid, p, N, M;

	int *A, *B, *x, *y, *z,*A_t, *B_t, *C_t;
	double time1, time2;
	//int Q;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &gpid);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	
	

	if(argc == 2)
	{
		if(gpid == 0) printf("\n\nnumP=%d, N=%s\n\n", p, argv[1]);
		M = atoi(argv[1]);
	}
	else
	{
		if(gpid == 0)
		{
			printf("Command Line : %s", argv[0]);
			for(i = 1; i < argc; i++)
			{
				printf("%s ", argv[i]);
			}
		}
		MPI_Finalize();
		return 0;
	}

	//M = N*p;
	N = M/p;
	//Q = M;

	A = (int*) malloc(N*M*sizeof(int));
	B = (int*) malloc(M*N*sizeof(int));
	x = (int*) malloc(M*M*sizeof(int));
	y = (int*) malloc(N*M*sizeof(int));
	B_t = (int*) malloc(M*N*sizeof(int));
	
	if(gpid == 0)	
	z = (int*) malloc(M*M*sizeof(int));

	

	initSubMatrix(A, N*M);// initialize row of A
	initSubMatrix(B, N*M);// initialize cols of B



	
	MPI_Barrier(MPI_COMM_WORLD);


	//Q = N*p;
	MPI_Request reqs[16];
	MPI_Status stats[16];
	int shift = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	time2 = MPI_Wtime();

	
        do
	{
	mulMatrix(A, B, y, N, M, shift,N, gpid); // local multiplication
	MPI_Irecv(B_t, (M*M/p), MPI_INT, (gpid+1)%p,(gpid+1)%p, MPI_COMM_WORLD, &reqs[gpid]);
	MPI_Send(B,   (M*M/p), MPI_INT, (gpid -1+p)%p, gpid, MPI_COMM_WORLD);

	MPI_Wait(&reqs[gpid], &stats[gpid]);
	*B = *B_t;
	shift++;
	
	
	
	} while(shift < p);
	MPI_Barrier(MPI_COMM_WORLD);
	time2 = MPI_Wtime() - time2;
	MPI_Gather(&y[0], N*M, MPI_INT, z, N*M, MPI_INT, 0, MPI_COMM_WORLD);// gather correctly the local product matrices


	

	//if(gpid == 0) printMatrix(z,N*p,M, "Z", gpid);
	
	if(gpid == 0) printf("%02d : Elapsed Time (PARALLEL) is %f\n", gpid, time2);
	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
