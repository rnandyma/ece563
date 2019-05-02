#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
struct node {
   float val;
   int node_num;
   struct node* l;
   struct node* r;
} nodeT;

#define MAXLEVEL 18

struct node* build(int level) {

   if (level < MAXLEVEL) {
      struct node* p = (struct node*) malloc(sizeof(nodeT));
      p->val = (float)rand()/RAND_MAX;
      p->node_num = level;
#pragma omp task 
{
      p->l = build(level+1);
}
#pragma omp task 
{
      p->r = build(level+1);
}
#pragma omp taskwait
      return p;
   } else {
      return NULL;
  }
}

/*void traverse(struct node* p, int* count_op) {
   //printf("the level is %d\n", p->node_num);
   if (p == NULL) return;
   if (p->l == NULL) return;
   else {
	if(p->node_num >=5 )
	{
	traverse(p->l, count_op);}
	else{
	#pragma omp task 
	{
	traverse(p->l, count_op);
	}	
	}
	}
   if (p->r == NULL) return;
   else {
	if(p->node_num >=5)
	{
	traverse(p->r, count_op);
	}
	else	
	{
	#pragma omp task 
        {	
	traverse(p->r, count_op);
	}	
	}
	}
#pragma omp taskwait
//#pragma omp critical
{
   if(p->val < 0.5)
   {
	(*count_op)++;
   }
 } 
return;
}*/

int traverse(struct node* p) {
   //printf("the level is %d\n", p->node_num);
   int count = 0;
   if (p == NULL) return count;
   if (p->l == NULL) return count ;
   else {
	if(p->node_num >=5 )
	{
	int count_tpex;
	count_tpex = traverse(p->l);
	//printf("the tpex is %d\n", count_tpex);
	 count = count+count_tpex;}
	else{
	#pragma omp task 
	{ int count_tpp;
	count_tpp = traverse(p->l);
	//printf("the tpp is %d\n", count_tpp);
	 count = count+count_tpp;
	}	
	}
	}
   if (p->r == NULL) return count;
   else {
	if(p->node_num >=5)
	{int count_tpex;
	count_tpex = traverse(p->r);
	//printf("the tpex is %d\n", count_tpex);
	count = count + count_tpex;
	}
	else	
	{
	#pragma omp task 
        {int count_tpp;
	count_tpp = traverse(p->r);
	//printf("the tpp is %d\n", count_tpp);	
	count = count + count_tpp;
	}	
	}
	}
#pragma omp taskwait
//#pragma omp critical
   if(p->val < 0.5)
   {
	count++;

   }
printf("the count is %d and the p_val is %f\n", count, p->val);
  
//#pragma omp taskwait
return count;
}



int main( ) {

double start, end;

#pragma omp parallel
{
#pragma omp single
{
   struct node* h = build(0);
  int count_op;
   //int init_count = 0;
   //int *count_ptr = &init_count;
start = omp_get_wtime();
   count_op = traverse(h);
end = omp_get_wtime();
   printf("the value of tid is %d count is %d\n",omp_get_thread_num(),count_op);
}
}

printf("the value of the exec time is %lf\n", end-start);


}
