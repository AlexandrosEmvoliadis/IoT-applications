#include <stdlib.h>
#include <stdio.h>
#include "cblas.h"
#include <time.h>
#include "cblas_f77.h"
#include <omp.h>
#include <float.h>
#include <mpi.h>
#include <stdbool.h>
#include <math.h>


//Project for PDS in Post Graduate Programme,ECE AUTH//
//**********USER NEED TO PROVIDE A TEXT DOCUMENT AND AN INTEGER*************//
//EACH LINE REPRESENTS THE COORDINATIONS FOR EACH POINT//
//EACH COORDINATE MUST BE SEPARATED WITH A SPACE CHARACTER AND AN END LINE CHARACTER TO CREATE A NEW POINT//
//THE ALGORITHM DOES NOT CATEGORIZE POINTS TO A CLUSTER,IT JUST FIND THE k NEAREST POINTS TO EACH POINT!!!//


typedef struct knnresult {
        int *nidx;
        double *ndist;
        int m;
        int k;
}knnresult;



void swap_id(int* a, int* b)
{
        int t = *a;
        *a = *b;
        *b = t;
}

void swap_dist(double* a, double* b)
{
        double t = *a;
        *a = *b;
        *b = t;
}

int partition (double *arr,int *id, int low, int high)
{
        double pivot = *(arr+low);
        int i = low;
        int j;

        for (j = low+1; j <= high; j++)
        {

                if (*(arr+j) <= pivot)
                {
                        i++;
                        swap_dist((arr+i),(arr+j));
                        swap_id((id+i),(id+j));
                }
        }
        swap_dist((arr+i),(arr+low));
        swap_id((id+i),(id+low));
        return i;
}

void quickSelect(double *arr,int *id, int low, int high, int k)
{
        if (low < high)
        {

                int pi = partition(arr,id,low, high);

                if(k <= pi) {
                        quickSelect(arr,id,low,pi - 1, k);
                }
                else if(k == pi+1) {
                        return;
                }
                else{
                        quickSelect(arr,id,pi + 1, high, k);
                }

        }
}

void merge(double* Dist_block,int* id_of_incoming,int number_of_points,int block_elements,int k_nearest,knnresult k)
{

        for(int i =0; i<number_of_points; i++) {
                for(int j = 0; j<block_elements; j++) {
                        if(isnan(*(Dist_block+j+i*block_elements))||*(Dist_block+j+i*block_elements)==INFINITY||*(Dist_block+j+i*block_elements)==-INFINITY) {
                                *(Dist_block+j+i*block_elements)=DBL_MAX;
                        }
                }
        }

        knnresult *merged = (knnresult*)malloc(sizeof(knnresult));
        merged->ndist = (double*)malloc(sizeof(double)*(number_of_points*(k_nearest+block_elements)));
        merged->nidx  = (int*)malloc(sizeof(int)*(number_of_points*(k_nearest+block_elements)));
        for (int i=0; i<number_of_points; i++) {
                for (int j=0; j<k_nearest+block_elements; j++) {
                        if(j<k_nearest) {
                                *(merged->ndist+j+(k_nearest+block_elements)*i) = *(k.ndist+j+k_nearest*i);
                                *(merged->nidx+j+(k_nearest+block_elements)*i) = *(k.nidx+j+k_nearest*i);
                        }
                        else{
                                *(merged->ndist+j+(k_nearest+block_elements)*i)=*(Dist_block+j-k_nearest+i*block_elements);
                                *(merged->nidx+j+(k_nearest+block_elements)*i) = *(id_of_incoming+j-k_nearest);
                        }
                }
                quickSelect((merged->ndist+(k_nearest+block_elements)*i),(merged->nidx+(k_nearest+block_elements)*i),0,k_nearest+block_elements-1,k_nearest);

                for (int j = 0; j<k_nearest; j++) {
                        *(k.ndist+j+k_nearest*i) = *(merged->ndist+j+(k_nearest+block_elements)*i);
                        *(k.nidx+j+k_nearest*i) = *(merged->nidx+j+(k_nearest+block_elements)*i);

                }


        }
        free(merged->nidx);
        free(merged->ndist);

}

void kNN(double *points1,double *points2,int* id1,int* id2,int number_of_points,int number_of_dimensions,int k_nearest, knnresult k1)
{
        int blocks = number_of_points/k_nearest;
        int remaining_elements = number_of_points%k_nearest;
        double *points1_squared = (double*)malloc(sizeof(double)*number_of_points);
        double *points2_squared = (double*)malloc(sizeof(double)*number_of_points);
        double *Dist = (double*)malloc(sizeof(double)*number_of_points*k_nearest);//distance matrix,blocked//
        for (int i=0; i<number_of_points; i++) {
                for (int j=0; j<number_of_dimensions; j++) {
                        *(points1_squared+i) += *(points1+j+i*number_of_dimensions) * *(points1+j+i*number_of_dimensions);
                        *(points2_squared+i) += *(points2+j+i*number_of_dimensions) * *(points2+j+i*number_of_dimensions);
                }
        }
        double alpha=-2.0;
        double beta=1.0;
        for (int b = 0; b<blocks; b++) {
                for (int i=0; i<number_of_points; i++) {
                        for (int j=0; j<k_nearest; j++) {
                                *(Dist+j+i*(k_nearest)) = *(points1_squared+i) + *(points2_squared+j+b*k_nearest);
                                if(*(id1+i)==*(id2+j+b*k_nearest)) {
                                        *(Dist+j+i*(k_nearest))=DBL_MAX;
                                }

                        }
                }
                //producing Distance_Matrix in blocks of k_nearest//
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,number_of_points,k_nearest,number_of_dimensions,alpha,points1,number_of_dimensions,(points2+b*k_nearest*number_of_dimensions),number_of_dimensions,beta,Dist,k_nearest);
                merge(Dist,id2+b*k_nearest,number_of_points,k_nearest,k_nearest,k1);
        }
        free(Dist);
        if (remaining_elements!=0) {
                if(blocks!=0) {
                        double *Dist_remaining =(double*)malloc(sizeof(double)*number_of_points*remaining_elements);
                        for (int i = 0; i<number_of_points; i++) {
                                for (int j=0; j<remaining_elements; j++) {
                                        *(Dist_remaining+j+i*remaining_elements) = *(points1_squared+i) + *(points2_squared+j+number_of_points-remaining_elements);
                                        if(*(id1+i)==*(id2+j+blocks*k_nearest)) {
                                                *(Dist_remaining+j+i*remaining_elements)=DBL_MAX;
                                        }
                                }
                        }
                        //producing Distance_Matrix with the remaining elements//
                        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,number_of_points,remaining_elements,number_of_dimensions,alpha,points1,number_of_dimensions,(points2+number_of_dimensions*(number_of_points-remaining_elements)),number_of_dimensions,beta,Dist_remaining,remaining_elements);
                        merge(Dist_remaining,(id2+blocks*k_nearest),number_of_points,remaining_elements,k_nearest,k1);
                        free(Dist_remaining);
                }
                else{
                        double *Dist_remaining =(double*)malloc(sizeof(double)*number_of_points*number_of_points);
                        for (int i = 0; i<number_of_points; i++) {
                                for (int j=0; j<number_of_points; j++) {
                                        *(Dist_remaining+j+i*number_of_points) = *(points1_squared+i) + *(points2_squared+j);
                                        if(*(id1+i)==*(id2+j)) {
                                                *(Dist_remaining+j+i*number_of_points)=DBL_MAX;
                                        }
                                }
                        }
                        //producing Distance_Matrix with the remaining elements//
                        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,number_of_points,number_of_points,number_of_dimensions,alpha,points1,number_of_dimensions,(points2),number_of_dimensions,beta,Dist_remaining,number_of_points);
                        merge(Dist_remaining,(id2),number_of_points,number_of_points,k_nearest,k1);
                        free(Dist_remaining);
                }
        }
        free(points1_squared);
        free(points2_squared);



}


int main(int argc, char* argv[])

{


        MPI_Init(NULL,NULL);
        int k = atoi(argv[2]);


        int rank,size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Status stat;
        MPI_Request req[2];
        int i=0;
        int ch=0;
        int lines=0;
        double X;
        int j=0;
        char chr;
        int rem;
        int ELEM[3];
        double *a,*QUERRY;


//GET ELEMENTS FROM TEXT IN ROOT PROCESS//
        if (rank == 0)
        {

                FILE *x = fopen(argv[1],"r");
                while(!feof(x))
                {
                        ch = fgetc(x);
                        if(ch == '\n')
                        {
                                lines++;
                        }
                        else if(ch==' ') {
                                i++;
                        }
                }
                fclose(x);
                x = fopen(argv[1],"r");


                int rows = i/lines;
                int rem = lines%size;

                a = (double*)malloc(lines*rows*sizeof(double));
                i = 0;
                j = 0;
                while(!feof(x)&&i<lines*rows) {

                        ch = fgetc(x);
                        if (ch!='\n'||ch!=' ') {
                                fscanf(x,"%lf",(a+i));
                                if(ch=='-') {
                                        *(a+i) = -*(a+i);
                                }


                        }
                        i++;



                }
                fclose(x);

                int ELEM[] = {lines,rows,rem};




        }

        //BROADCAST DIMENSIONS TO ALL//
        MPI_Bcast(&ELEM,3,MPI_INT,0,MPI_COMM_WORLD);
        rem = ELEM[2];
        int D = ELEM[1];
        int s = ELEM[0]/size;
        if (rank!=0) {
                a = (double*)malloc(ELEM[0]*ELEM[1]*sizeof(double));

        }
        //BROADCAST TEXT TO ALL//
        MPI_Bcast(a,ELEM[0]*ELEM[1],MPI_DOUBLE,0,MPI_COMM_WORLD);
        int *loc = malloc(sizeof(int)*size);
        int*disp =  malloc(sizeof(int)*size);
        int elem_check = 0;
        if (rem!=0) {
                s++;
                if(rank>=rem) {
                        elem_check = 1;
                }

        }
        s = D*s;

        int N = s/D;

        int prev,next;

        if (rank==0) {
                prev = size-1;
                next = rank+1;
        }
        else if(rank==size-1) {
                next = 0;
                prev = rank-1;
        }
        else{
                next = rank+1;
                prev = rank-1;
        }

        int K;
        if(elem_check==1) {
                K = N-1;
        }
        else{
                K=N;
        }
        K= K*D;

        MPI_Gather(&K,1,MPI_INT,loc,1,MPI_INT,0,MPI_COMM_WORLD);

        if(rank==0) {
                int sum = 0;
                for (i=0; i<size; i++) {
                        disp[i] = sum;
                        sum+=loc[i];


                }
        }

        MPI_Bcast(disp,size,MPI_INT,0,MPI_COMM_WORLD);
        QUERRY = (double*)malloc(sizeof(double)*N*D);
        int *MY_ID = (int*)malloc(sizeof(int)*N);

        for(i=0; i<N; i++) {
                for(int j=0; j<D; j++) {
                        *(QUERRY+i*D+j) = *(a+i*D+j+disp[rank]);

                }

                *(MY_ID+i) = (i+rank*N);


        }
        free(a);
        free(loc);
        free(disp);
        if(elem_check==1) {
                for(i=0; i<D; i++) {
                        *(QUERRY+s-1-i) = DBL_MAX;
                }
        }
        /*initialize k1*/
        knnresult *k1 = (knnresult*)malloc(sizeof(knnresult));
        k1->nidx = (int*)malloc(sizeof(int)*N*k);
        k1->ndist = (double*)malloc(sizeof(double)*N*k);
        for (int i = 0; i<N; i++) {
                for (int j = 0; j<k; j++) {
                        *(k1->nidx+k*i+j) = -1;
                        *(k1->ndist+k*i+j) = DBL_MAX;
                }
        }

        double* c1 = (double*)malloc(sizeof(double)*N*D);
        int* id   = (int*)malloc(sizeof(int)*N);

        c1= QUERRY;
        id = MY_ID;
        double start,end;

        if (rank==0) {
                start = MPI_Wtime();
        }
        MPI_Isend(c1,N*D,MPI_DOUBLE,next,rank,MPI_COMM_WORLD, &req[0]);
        MPI_Irecv(c1,N*D,MPI_DOUBLE,prev,prev,MPI_COMM_WORLD, &req[1]);
        MPI_Isend(id,N,MPI_INT,next,rank,MPI_COMM_WORLD, &req[0]);
        MPI_Irecv(id,N,MPI_INT,prev,prev,MPI_COMM_WORLD, &req[1]);
        kNN(QUERRY,c1,MY_ID,id,N,D,k,*k1);
        MPI_Waitall(2,req,MPI_STATUSES_IGNORE);

        for (int itt=0; itt<size-1; itt++)
        {

                MPI_Waitall(2,req,MPI_STATUSES_IGNORE);

                MPI_Isend(c1,N*D,MPI_DOUBLE,next,rank,MPI_COMM_WORLD,&req[0]);
                MPI_Irecv(c1,N*D,MPI_DOUBLE,prev,prev,MPI_COMM_WORLD,&req[1]);
                MPI_Isend(id,N,MPI_INT,next,rank,MPI_COMM_WORLD, &req[0]);
                MPI_Irecv(id,N,MPI_INT,prev,prev,MPI_COMM_WORLD, &req[1]);


                kNN(QUERRY,c1,MY_ID,id,N,D,k,*k1);
                if (rank==0) {
                        printf("%d\n",itt );
                }

        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank==0) {
                end = MPI_Wtime();
                printf("Time elapsed during the job: %.2fs.\n", end - start);
        }






        MPI_Finalize();

        return 0;
}
