/****************************************************************************
 ****************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
extern void draw_heat(int nx, int ny);       /* X routine to create graph */

#define M      360                       /* x dimension of problem grid /Collumns*/
#define N      360                       /* y dimension of problem grid /Rows */
#define ITER       300000000              /* number of time steps */
#define MASTER      0                  /* taskid of first process */
#define TAM_BLOCO 122

#define TOL tol
#define CLUSTER 100
#define TIME_RESOLUTION 1000000 /* time measuring resolution (us) */
#define epsilon FLT_EPSILON




double u[N][N], w[N][N];
double tol = 1.0f / (N * N);

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};








/*****************************************************************************
 *  subroutine init
 *****************************************************************************/

void init(){
    int i, j;
    for (j = 0; j < N; j++) u[N-1][j] = 100.0f; //bottom
    for (j = 0; j < N; j++) u[0][j] = 0.0f; //top
    for (i = 0; i < N; i++) u[i][N-1] = 100.0f; //right
    for (i = 0; i < N; i++) u[i][0] = 100.0f;//left

    for (j = 0; j < N; j++) w[N-1][j] = 100.0f; //bottom
    for (j = 0; j < N; j++) w[0][j] = 0.0f; //top
    for (i = 0; i < N; i++) w[i][N-1] = 100.0f; //right
    for (i = 0; i < N; i++) w[i][0] = 100.0f;//left

    for (i = 1; i < (N-1); i++)
        for (j = 1; j < (N-1); j++) u[i][j] = w[i][j] = 50.0f;

}



/**************************************************************************
 * subroutine compute por blocos de linhas
 **************************************************************************/


void compute(double* linhas, double * linha, int linesByBlock){
  int i,j;
  //#pragma omp parallel for private(j)
  for (i=1;i<linesByBlock-1;i++){
    for (i = 1; i < linesByBlock-1; i++)
                for (j = (i %2 == 0? 2 : 1); j < (N-1); j+=2){
      linhas[j+((i-1)*M)]= 0.25*(linhas[((i-1)*M)+j]+linhas[((i+1)*M)+j]/*calcular a colum */
                      + linhas[(i*M)+j-1] + linhas[(i*M)+j+1] /*calulate line*/       
                );
     linha[j+((i-1)*M)] = linhas[j+((i-1)*M)];
    }
    //#pragma omp parallel for private(j)
     for (i = 1; i < linesByBlock-1; i++){
     linha[((i-1)*M)]=100;//heat font
            for (j = (i %2 == 0? 1 : 2); j < (N-1); j+=2){
      linha[j+((i-1)*M)]= 0.25*(linhas[((i-1)*M)+j]+linhas[((i+1)*M)+j]/*calcular a colum */
                      + linhas[(i*M)+j-1] + linhas[(i*M)+j+1] /*calulate line*/       
                );
    }
  linha[j+((i-1)*M)]=100;//heat font
}
}
}


/**************************************************************************
 * WORK
 **************************************************************************/

/*
int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)


*/


int main (int argc, char *argv[]){

int taskid,                     /* this task's unique id */
  numtasks,                   /* number of tasks */
  source,               /* to - from for message send-receive */
  index,              /* loop variables */
    proc;
int iter;
MPI_Status status;
int j,i;
double ** temp;
   int blocotam = TAM_BLOCO;
double diff = TOL +1, auxDiff;
double linha[M*(blocotam-2)];//, enviada pelos slaves;

double enviada [blocotam*M];
int aux, nrslave;
double time;

/* First, find out my taskid and how many tasks are running */
   MPI_Init(NULL,NULL);
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
   time =  MPI_Wtime(); 



   //printf("taskid:%d\n",taskid);
   if (taskid == MASTER) {
      /************************* master code *******************************/
      /* Initialize & fill matrix */
      init();

      for (iter = 0; iter < ITER; iter++)
      {
          printf("====Iteração nr %d===\n", iter);

          #pragma omp parallel for private(i, aux, j)
          for (index = 1, nrslave = 1; index < N; index = index + (blocotam - 2), nrslave++)
          {

              for (i = index - 1, aux = 0; i < index + blocotam - 2 && i < N; i++, aux++)
              { //CONFIRMAR ESTE IF

                  for (j = 0; j < M; j++)
                  {
                      enviada[(aux * M) + j] = w[i][j];
                      //printf("Valor enviado w[%d][%d]=%f com index %d \n", i,j,w[i][j],index);
                  }
              }

              MPI_Send(enviada, blocotam * M, MPI_DOUBLE, nrslave, index, MPI_COMM_WORLD);
          }
          /* Now wait for results from all worker tasks */

          // tag vai conter a linha onde o Master vai inserir em u
          int x = N / (blocotam - 2);
          for (; x > 0; x--)
          {

              MPI_Recv(linha, (M * (blocotam - 2)), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

              //Atualization of u
              int line;
              line = status.MPI_TAG;
              for (i = 0; i < (blocotam - 2) * M && line < N - 1; line++)
              {
                  for (j = 0; j < M; j++, i++)
                  {    
                    auxDiff = (linha[i] - w[line][j]);
                    diff = diff > aux ?  : auxDiff;
                    printf(" diff: %lf \n",diff);
                    w[line][j] = linha[i];
                  }
              }
          }
         if ( diff <= TOL) {
             printf ("Iteração : %d", iter);
            time =  MPI_Wtime() - time;
            printf("Elapsed time %lf\n",time);
             MPI_Finalize();
         }
      }
   }  
 
   /* End of master code */

   /************************* workers code **********************************/
  


  else{

    for (i=0;i<ITER;i++){

      MPI_Recv(enviada, (M*blocotam), MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 

        compute(enviada, linha, blocotam);

        MPI_Send(linha,M*(blocotam-2), MPI_DOUBLE,0,status.MPI_TAG,MPI_COMM_WORLD); //check
        //printf("B\n");
        for(j=0;j<N*(blocotam-3);j++){
            //printf("Nodo : %d DEPOIS DO COMPUTE:linha[%d]=%lf\n", taskid,j,linha[j]);
        }
    
    }
  }
      MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
      MPI_Finalize();
}


