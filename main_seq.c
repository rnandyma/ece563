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
      //printf("the level is %d\n", p->node_num);
      //p->val = 1;

      p->l = build(level+1);
      p->r = build(level+1);

      return p;
   } else {
      return NULL;
  }
}

void traverse(struct node* p, int* count_op) {
   //printf("the level is %d\n", p->node_num);
   if (p == NULL) return;
   if (p->l == NULL) return;
   else {
	traverse(p->l, count_op);
	}
   if (p->r == NULL) return;
   else {	
	traverse(p->r, count_op);
	}
   if(p->val < 0.5)
   {
	(*count_op)++;
   }

   
return;
}

int main( ) {

double start, end;





   struct node* h = build(0);
  // int count_op;
   int init_count = 0;
   int *count_ptr = &init_count;
start = omp_get_wtime();
   /*count_op = */traverse(h,count_ptr);
end = omp_get_wtime();
   printf("the value of count is %d\n",*count_ptr);



printf("the value of the exec time is %lf\n", end-start);


}
