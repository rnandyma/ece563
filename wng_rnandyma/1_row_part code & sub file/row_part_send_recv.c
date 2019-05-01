#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void mulMatrix(int *A, int *B, int *C, int N, int M, int circular_shift, int block_width, int pid)
{
	int i, j, k;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			
			C[(i*M)+(j+(shift*block_width)+(pid*block_width))%M] = 0;
			for(k = 0; k < M; k++)
			{
				
				C[(i*M)+(j+(shift*block_width)+(pid*block_width))%M] += A[(M*i+k)] * B[(j*M+k)];
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
	double time1;
	


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

	
	N = M/p;
	

	A = (int*) malloc(N*M*sizeof(int)); // create matrix for row partioned matrix A
	B = (int*) malloc(M*N*sizeof(int)); // create matrix for column partitioned matrix B
	x = (int*) malloc(M*M*sizeof(int)); // final total matrix of size M*M
	y = (int*) malloc(N*M*sizeof(int)); // temporary product matrix in every process	
	B_t = (int*) malloc(M*N*sizeof(int)); // temporary column matrix in every process to receive the bunch of columns from other process
	
	if(gpid == 0)	
	z = (int*) malloc(M*M*sizeof(int));

	initSubMatrix(A, N*M);// initialize row of A
	initSubMatrix(B, N*M);// initialize cols of B
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Request reqs[16];
	MPI_Status stats[16];
	int circular_shift = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = MPI_Wtime();

	
    do
	{
	mulMatrix(A, B, y, N, M, circular_shift,N, gpid); // local multiplication
	MPI_Irecv(B_t, (M*M/p), MPI_INT, (gpid+1)%p,(gpid+1)%p, MPI_COMM_WORLD, &reqs[gpid]);
	MPI_Send(B,   (M*M/p), MPI_INT, (gpid -1+p)%p, gpid, MPI_COMM_WORLD);

	MPI_Wait(&reqs[gpid], &stats[gpid]);
	*B = *B_t;
	circular_shift++;
	
	
	
	} while(shift < p);
	MPI_Barrier(MPI_COMM_WORLD);
	time1 = MPI_Wtime() - time1;
	MPI_Gather(&y[0], N*M, MPI_INT, z, N*M, MPI_INT, 0, MPI_COMM_WORLD);// gather correctly the local product matrices


	

	//if(gpid == 0) printMatrix(z,N*p,M, "Z", gpid);
	
	if(gpid == 0) printf("%02d : Elapsed Time (PARALLEL) is %f\n", gpid, time2);
	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
