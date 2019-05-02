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

struct node* build_p(int level) {

   if (level < MAXLEVEL) {
      struct node* p = (struct node*) malloc(sizeof(nodeT));
      p->val = (float)rand()/RAND_MAX;
      p->node_num = level;
if(level < 5)
{
#pragma omp task 
{
      p->l = build_p(level+1);
}
#pragma omp task 
{
      p->r = build_p(level+1);
}
#pragma omp taskwait
      return p;
}

else
{
      p->l = build_p(level+1);
      p->r = build_p(level+1);
      return p;
}

   } else {
      return NULL;
  }
}

struct node* build_s(int level) {

   if (level < MAXLEVEL) {
      struct node* p = (struct node*) malloc(sizeof(nodeT));
      p->val = (float)rand()/RAND_MAX;
      p->node_num = level;
      //printf("the level is %d\n", p->node_num);
      //p->val = 1;

      p->l = build_s(level+1);
      p->r = build_s(level+1);

      return p;
   } else {
      return NULL;
  }
}



int traverse_p(struct node* p) {
   //printf("the level is %d\n", p->node_num);
   int count_left = 0;
   int count_right = 0;
   if (p == NULL) return (count_left+count_right);
   if (p->l == NULL) return (count_left+count_right) ;
   else {
	if(p->node_num < 5){
	#pragma omp task shared(count_left)
	count_left = traverse_p(p->l);}
	else {
	count_left = traverse_p(p->l);	
	}
	
	}
   if (p->r == NULL) return (count_left+count_right);
   else {
	if(p->node_num < 5){
	#pragma omp task shared(count_right)
	count_right = traverse_p(p->r);}
	else{
	count_right = traverse_p(p->r);
	}
	}
if(p->node_num < 5)
{
#pragma omp taskwait
//#pragma omp critical
   if(p->val < 0.5)
   {
	return (count_left+count_right+1);

   }
  else
	return (count_left+count_right);
}
else
{
if(p->val < 0.5)
   {
	return (count_left+count_right+1);

   }
  else
	return (count_left+count_right);
}

}

int traverse_s(struct node* p) {
   //printf("the level is %d\n", p->node_num);
   int count_left = 0;
   int count_right = 0;
   if (p == NULL) return (count_left+count_right);
   if (p->l == NULL) return (count_left+count_right);
   else {
	count_left = traverse_s(p->l);
	}
   if (p->r == NULL) return (count_left+count_right);
   else {	
	count_right = traverse_s(p->r);
	}
   if(p->val < 0.5)
   {
	return (count_left+count_right+1);
   }
   else
	return count_left + count_right;
}



int main( ) {

double start_s, end_s;
double start_p, end_p;

#pragma omp parallel
{
#pragma omp single
{
    int count_opp, count_ops;
    
start_s = omp_get_wtime();   
  struct node* h = build_p(0);
  
   count_ops = traverse_p(h);
end_s = omp_get_wtime();
   //printf("the value of tid is %d count is %d\n",omp_get_thread_num(),count_ops);
   printf("the value of the exec time seq is %lf\n", end_s-start_s);
start_p = omp_get_wtime();
   count_opp = traverse_p(h);
end_p = omp_get_wtime();
   //printf("the value of tid is %d count is %d\n",omp_get_thread_num(),count_opp);
   //printf("the value of the exec time parallel is %lf\n", end_p-start_p);
}
}




}
