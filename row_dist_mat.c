#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
void mulMatrix(int *A, int *B, int *C, int N, int M, int Q)
{
	int i, j, k;

	for(i = 0; i < N; i++)
	{
		for(j = 0; j < Q; j++)
		{
			C[(Q*i+j)] = 0;
			for(k = 0; k < M; k++)
			{
				C[(Q*i+j)] += A[(M*i+k)] * B[(Q*k+j)];
			}
		}
	}

	return;
}

void initVector(int *A, int N)
{
	int i;

	for(i = 0; i < N; i++)
	{
		A[i] = rand()%3;
	}
	return;
}

void printVector(int *A, int N, int M, char *str, int gpid)
{
	int i, j;

	printf("%s (%02d)\n", str, gpid);
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < M; j++) printf("%d ", A[(M*i+j)]);
	   printf("\n");
	}

	return;
}

int main(int argc, char* argv[])
{
	int i, j, k;
	int pid, nump, N, M;

	int *A, *B, *x, *y, *A_t, *B_t, *C_t;
	double time1, time2;
	int Q;


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &nump);

	MPI_Datatype sblock, sblocktype;
	int disp[p]; // to be used in gatherv and scatterv, might remove it if I use only MPI_Send and recv
	int scount[p]; // to be used in gatherv and scatterv, might remove it if I use only MPI_Send and recv
	int rcount[p]; // to be used in gatherv and scatterv, might remove it if I use only MPI_Send and recv
	
	srand();

	if(argc == 2)
	{
		if(gpid == 0) printf("\n\nnumP=%d, N=%s\n\n", p, argv[1]);
		N = atoi(argv[1]);
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
