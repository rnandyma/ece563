//#define N 5
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stddef.h>
#include "mpi.h"


//void print_results(char *prompt, int a[N][N]);

int main(int argc, char *argv[])
{
    int i, j, k, pid, numP, sum = 0, rc = 0;
    double time1, time2;
    int l ;
    int N;
    N = atoi(argv[1]);
    int a[N][N];
    int b[N][N];
    //int a[N][N]={{1,2,3,4,4},{5,6,7,8,8},{9,1,2,3,3},{4,5,6,7,7},{8,9,10,11,12,}};
    //int b[N][N]={{1,2,3,4,4},{5,6,7,8,8},{9,1,2,3,3},{4,5,6,7,7},{8,9,10,11,12,}};
    int output[N][N];
    int c[N][N],d[N][N];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numP);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Status stat;
    for (int x = 0; x < N; x++)
            {
                    for (int y = 0; y < N; y++)
                    {
                         //output[x][y] = 0;
			  a[x][y] = 1;
			  b[x][y] = 1;
			  c[x][y] = 0;
			  d[x][y] = 0;              
                    }
                   
            }


    

    //scatter rows of first matrix to different processes  
    for(l = 0; l<= N/numP; l++){   
    MPI_Scatter(a[l*numP], N, MPI_INT, c[l], N, MPI_INT,0,MPI_COMM_WORLD);}
    //scattering the rows of B or we can send columns too, better if columns
    /*for(l = 0; l<= N/numP; l++){   
    MPI_Scatter(a[l*numP], N, MPI_INT, c[l], N, MPI_INT,0,MPI_COMM_WORLD);}*/


    //broadcast second matrix to all processes
    //MPI_Bcast(b, N*N, MPI_INT, 0, MPI_COMM_WORLD);
    

    //MPI_Barrier(MPI_COMM_WORLD);
           for(int i = 0; i<=N/numP; i++){
          //perform vector multiplication by all processes
          for (j = 0; j < N; j++)
            {       output[i][j] = 0;
		printf("the value of c(%d %d) is %d\n", i,j,c[i][j]);
                    for (k = 0; k < N; k++)
                    {       
                            output[i][j] = output[i][j] + (c[i][k] * b[k][j]);
                                       
                    } //printf("the value of output(%d,%d) is %d\n", i,j,output[i][j]); 
                    //cc[k][i] = sum;
                    //sum = 0;
            }
MPI_Gather(output[i], N, MPI_INT, d[i*numP], N, MPI_INT, 0, MPI_COMM_WORLD);
}
    //for(int m = 0; m*size*N<N*N; m= m+1){
    //MPI_Gather((cc+m*N), N, MPI_INT, (c+m*size*N), N, MPI_INT, 0, MPI_COMM_WORLD);}

    MPI_Barrier(MPI_COMM_WORLD); 
    if (pid == 0)  {
     for(i = 0; i<N; i++){
	for(j = 0; j < N; j++)
        {printf("%d ", d[i][j]);}
	printf("\n");}

	}                     
    MPI_Finalize();
             //I_ADDED_THIS
        
}

