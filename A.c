#include  <stdio.h>
#include  <omp.h>

int main(int argc, char *argv[]) { 
 
   //char* n = argv[1];
   //int numThreads = 0;
   int numProcs ;

  /* for (int i = 0; n[i] != 0; i++) {
      numThreads = numThreads*10 + n[i] - '0';
   }*/

   //omp_set_num_threads(numThreads);
    numProcs = omp_get_num_procs();
    printf("Number of processors are %d\n",numProcs);

#pragma omp parallel
{
   printf("Hello world from thread %d\n",omp_get_thread_num( ));
   #pragma omp single
    printf("The thread executing the single statement is %d\n", omp_get_thread_num());
   #pragma omp master
    printf("The thread executing the master statement is %d\n", omp_get_thread_num());
}

   return 0;
}


