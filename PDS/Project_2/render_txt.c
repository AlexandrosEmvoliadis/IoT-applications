#include<stdlib.h>
#include<stdio.h>
#include"cblas.h"
#include<time.h>
#include "cblas_f77.h"
#include <omp.h>
#include <float.h>

double *renderX(int a,int d){

 double upper = 100.0;
 double lower = -100.0;
 double *A = (double*)malloc(a*d*sizeof(double));
 srand(time(NULL));
 for (int i = 0;i<a;i++){
    for(int j = 0;j<d;j++){
        *(A+d*i+j) =  (double)rand()/(RAND_MAX/(double)(upper-lower))+lower;
    }
 }

 return A;
}
int main(){
  int N,D;
double *X;
printf("Give the Number of Points in X");
    scanf("%d", &N);
    printf("Give the Number of Dimensions ");
    scanf("%d", &D);

X = renderX(N,D);
FILE *x;
x = fopen("X_random.txt","w");
for (int i = 0;i<N;i++){
  for (int j = 0;j<D;j++){

    fprintf(x,"%lf ",*(X+D*i+j));
  }
  fprintf(x,"\n");
}
fclose(x);
free(X);
return 0;
}
