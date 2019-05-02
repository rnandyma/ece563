#include <stdio.h>
#include <stdlib.h>
#include <math.h>
struct proc_info
{
int low_index;
int high_index;
int difference;
}nodeP;
int main(int argc, char *argv[])
{

printf("//////////////////////////////////////////////new data input/////////////////////////////////\n");
int numP, array_size, chunk_size;
numP = atoi(argv[1]);
array_size = atoi(argv[2]);
chunk_size = atoi(argv[3]);

struct proc_info proc_array[numP];

int array[array_size];




for(int i = 0; i<array_size; i++)
array[i]=1;

printf("the number of procs is %d , array of size %d, chunk_size for block_cyclic is %d\n", numP, array_size, chunk_size);
///////////////////////////for cyclic distribution//////////////////////////////////////////
printf("/////////////////////////////////////////////////cyclic distribution/////////////////////////////////////////////\n");
for(int l = 0; l<numP; l++)
proc_array[l].low_index = l;

int k = 0;
for (int j = numP; j<array_size; j++)
{
proc_array[k].high_index = j;
proc_array[k].difference = numP;

k++;
if(k == numP)
k = 0;
}
for(int m = 0; m<numP; m++)
printf("%d : [%d : %d : %d]\n", m,proc_array[m].low_index, proc_array[m].high_index, numP);



////////////////////////////for block distribution///////////////////////////////////////////////////
printf("/////////////////////////////////////////block distribution//////////////////////////////////////\n");
int remainder = array_size%numP;
for(int a = 0; a<numP; a++)
{
if(a == 0)
proc_array[a].low_index = 0;
else 
proc_array[a].low_index = proc_array[a-1].high_index + 1;
if(remainder == 0)// all blocks will be of same size
proc_array[a].high_index = proc_array[a].low_index + (array_size/numP) -1;
else // first remainder blocks will have ceil(n/p) elements rest numP-remaninder will have floor(n/p) elements
{
if(a <remainder)
proc_array[a].high_index = proc_array[a].low_index + (array_size/numP);
else
proc_array[a].high_index = proc_array[a].low_index + (array_size/numP) -1;
}
/*else{
if((a%2 == 0) && (a!=0))
proc_array[a].high_index = proc_array[a].low_index + (array_size/numP);
else 
proc_array[a].high_index = proc_array[a].low_index + (array_size/numP) -1;}*/
proc_array[a].difference = 1;
}
for(int m = 0; m<numP; m++)
printf("%d : [%d : %d : %d]\n", m,proc_array[m].low_index, proc_array[m].high_index, proc_array[m].difference);


//////////////////////////////for block cyclic distribution /////////////////////////////////////////////////
printf("/////////////////////////////////////////block cyclic distribution//////////////////////////////////////////\n");
int row_size = (ceil((ceil(array_size/chunk_size))/numP))*chunk_size;
//printf("row_size = %d\n", row_size);

int array_proc[numP][row_size];
for(int i =0;i<numP;i++)
{
for(int j=0;j<row_size;j++)
array_proc[i][j] = -1;
}


int b = 0;
int count = 0;
for( int n = 0; n <array_size; n=n+chunk_size)
{
for(int z = n;z<n+chunk_size; z++)
{array_proc[b][(z-n)+(count*chunk_size)] = z;/*printf("processor %d at col %d contains the element of index %d\n",b,z,array_proc[b][z]);*/}

b++;
if(b==numP)
{b = 0; count++;}

}
for(int i =0;i<numP;i++)
{
for(int j=0;j<row_size;j++)
{
if(array_proc[i][j] != -1)
printf("processor %d contains the element of index %d\n",i,array_proc[i][j]);
}
}
}
