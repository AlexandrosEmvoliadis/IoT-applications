#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "mmio.c"
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <cilk/cilk.h>



uint32_t M, N,n, nz, i, j, k,l, *coo_col, *coo_row,sum ,trngls;




int main(int argc, char *argv[]) {



 int ret_code;
 MM_typecode matcode;
 FILE *f;
 double *val;

 if (argc < 2)
{
 fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
 exit(1);
}
 else
 {
     if ((f = fopen(argv[1], "r")) == NULL)
         exit(1);
 }

 if (mm_read_banner(f, &matcode) != 0)
 {
     printf("Could not process Matrix Market banner.\n");
     exit(1);
 }

 if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
         mm_is_sparse(matcode) )
 {
     printf("Sorry, this application does not support ");
     printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
     exit(1);
 }

 if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
     exit(1);


 coo_col = (uint32_t *) malloc(nz * sizeof(uint32_t));
 coo_row = (uint32_t *) malloc(nz * sizeof(uint32_t));
 val = (double *) malloc(nz * sizeof(double));


 if (!mm_is_pattern(matcode))
 {
 for (i=0; i<nz; i++)
 {
     fscanf(f, "%d %d %lg\n", &coo_row[i], &coo_col[i], &val[i]);
     coo_row[i]--;
     coo_col[i]--;
 }
 }
 else
 {
 for (i=0; i<nz; i++)
 {
     fscanf(f, "%d %d\n", &coo_row[i], &coo_col[i]);
     val[i]=1;
     coo_row[i]--;
     coo_col[i]--;
 }
 }

 if (f !=stdin) fclose(f);
 mm_write_banner(stdout, matcode);
 mm_write_mtx_crd_size(stdout, M, N, nz);

 uint32_t* col_pop  = (uint32_t*)malloc((M+1)*sizeof(uint32_t));
 for (i = 0;i<M+1;i++){
      col_pop[i] = 0;
    }

    /*coo2csr happens here*/
   uint32_t sum  = 0;
   int i = 0;
    while (i<M+1){
      col_pop[i] = sum;
      while (coo_col[sum]==i){

            sum++;
        }

        i++;

    }
    col_pop[M] = nz;
    free(coo_col);/*useless*/

    printf("\n");
    int32_t * c3 = (uint32_t *)malloc(M * sizeof(uint32_t));
    /*initialize the number of triangles each node participates in*/
     for (i=0;i<M;i++){
       c3[i]=0;
     }

     double time_spent = 0.00000;
     trngls = 0;
     clock_t begin = clock();
     /*counting triangles*/
     for(i=0;i<M;i++){
       if (col_pop[i+1]-col_pop[i]>1){/*DO NOT GET INSIDE A POSITION CONNECTED WITH NO MORE THAN 1 ELEMENTS*/
        for(j=col_pop[i];j<col_pop[i+1];j++){
          for(k=col_pop[coo_row[j]];k<col_pop[coo_row[j]+1];k++){
            for(l=j+1;l<col_pop[i+1];l++){
              if ( coo_row[l]==coo_row[k]){
                c3[i]++;
                c3[coo_row[l]]++;
                c3[coo_row[j]]++;
                trngls = trngls+3;

              }
            }
          }
        }
      }
    }
    clock_t finish = clock();
    time_spent += (double)(finish - begin) / CLOCKS_PER_SEC;
    free(col_pop);
    free(coo_row);
    free(c3);

    printf("The elapsed time is %f seconds", time_spent);

    printf("%d\n",trngls );


return 0;

}
