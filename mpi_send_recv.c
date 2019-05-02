/*  C  */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

// For Block (Second method) distribution
#define BSIZE(pid, numP, N)           (((pid+1)*N/numP) - (pid*N/numP))
#define BINIT(pid, numP, N, i)        ((pid*N/numP) + i)

int main(int argc, char* argv[])
{
   int pid;
   int numP, N;
   int i, j, alen, t, sum, rc;

   int* A;

   double time1, time2;

   MPI_Request reqs[16];
   MPI_Status stats[16];


   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numP); 
   MPI_Comm_rank(MPI_COMM_WORLD, &pid);

   N = atoi(argv[1]);

   sum = 0;
   MPI_Barrier(MPI_COMM_WORLD);

   
      alen = BSIZE(pid, numP, N);
      A = (int*) malloc(alen*sizeof(int));
      for(i = 0; i < alen; i++)
      {
         A[i] = 1;
      }
   

   
   MPI_Barrier(MPI_COMM_WORLD);

   
   time2 = MPI_Wtime();
   MPI_Barrier(MPI_COMM_WORLD);

   
      for(j = 0; j < 1000; j++)
      {
         t = 0;
         for(i = 0; i < alen; i++) t += A[i];
       
         if(pid%2 == 1) rc = MPI_Send(&t, 1, MPI_INT, (pid-1), pid, MPI_COMM_WORLD);
         else
         {
            rc = MPI_Recv(&sum, 1, MPI_INT, (pid+1), (pid+1), MPI_COMM_WORLD, &stats[pid]);
            //if(j == 0) //printf("Recv %02d : %d\n", pid, t);
            t += sum;
            //if(j == 0) //printf("Recv %02d : %d\n", pid, t);
         }
         
         if(pid%4 == 2) rc = MPI_Send(&t, 1, MPI_INT, (pid - 2), pid, MPI_COMM_WORLD);
         else if(pid%4 == 0)
         {
            rc = MPI_Recv(&sum, 1, MPI_INT, (pid+2), (pid+2), MPI_COMM_WORLD, &stats[pid]);
            t += sum;
            //if(j == 0) //printf("Recv %02d : %d\n", pid, t);
         }

         if(pid%8 == 4) rc = MPI_Send(&t, 1, MPI_INT, (pid-4), pid, MPI_COMM_WORLD);
         else if(pid%8 == 0)
         {
            rc = MPI_Recv(&sum, 1, MPI_INT, (pid+4), (pid+4), MPI_COMM_WORLD, &stats[pid]);
            t += sum;
            //if(j == 0) //printf("Recv %02d : %d\n", pid, t);
         }

         if(pid%16 == 8) rc = MPI_Send(&t, 1, MPI_INT, (pid-8), pid, MPI_COMM_WORLD);
         else if(pid%16 == 0)
         {
            rc = MPI_Recv(&sum, 1, MPI_INT, (pid+8), (pid+8), MPI_COMM_WORLD, &stats[pid]);
            t += sum;
            //if(j == 0) //printf("Recv %02d : %d\n", pid, t);
         }
      }
   

   MPI_Barrier(MPI_COMM_WORLD);   
   time1 = (MPI_Wtime() - time2)/1000;

   if(pid == 0) printf("(Send/Recv)  Time Elapsed (averaged over 1000 runs) for size %d  : %f\n", N, time1);
   if(pid == 0) printf("pid : %02d,   sum : %03d\n", pid, t);

   MPI_Barrier(MPI_COMM_WORLD);   
   MPI_Finalize();

   return 0;
}
