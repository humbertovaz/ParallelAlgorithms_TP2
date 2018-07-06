#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
extern void draw_heat(int nx, int ny);       /* X routine to create graph */

#define M      200                       /* x dimension of problem grid /Collumns*/
#define N      200                       /* y dimension of problem grid /Rows */
#define ITER       3                 /* number of time steps */
#define NSLAVES 3
#define MASTER      0                  /* taskid of first process */

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
 * subroutine compute
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
double diff = TOL +1, auxDiff;
double linha[M];//, enviada;
int taskid,                     /* this task's unique id */
  numtasks,                   /* number of tasks */
  source,               /* to - from for message send-receive */
  index,              /* loop variables */
    proc;
int iter;
MPI_Status status;
int j,i;
double enviada [3*M];
int aux, auxindex, env;
int slaves= NSLAVES;
int acabado;
double time;
double temp;
time = MPI_Wtime();
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
        /*  Now send startup information to each worker  */
        /* Need to send a contignuos struct */
        //double enviada [3*M];
      /* Distribute work to workers.  1 row per worker */
    for(iter=0; diff > TOL ;iter++){
        diff=0;
        //


 //PARTE ESTATICA
      for (env=1, index=1; index<N-1 && index <=slaves; env++, index++){

        //CREATE THE ARRAY THAT WILL BE SEND
        for(i=index-1, aux=0; i<index+2; i++,aux++){
          for (j=0;j<M;j++){
            enviada[(aux*M)+j]=w[i][j];
          }
        }

        acabado=0;
        MPI_Send(&acabado, 1, MPI_INT,env,index,MPI_COMM_WORLD);
        //
        MPI_Send(enviada, 3*M, MPI_DOUBLE,env,index,MPI_COMM_WORLD); //Envia um bloco de linhas para o iésimo slave 
      //
      }

      env--;
      index--;


  //PARTE DINAMICA
      for (auxindex=1; auxindex<N-1; auxindex++){
        // tag vai conter a linha onde o Master vai inserir em u 
        MPI_Recv(linha, M, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        //Atualization of u
        for(aux=0; aux < M;aux++){
            u[(status.MPI_TAG)][aux]= linha[aux];
        //    printf("Iteracao nr: %d u[%d][%d]=%lf\n",iter,status.MPI_TAG,i,u[status.MPI_TAG][i] );
        }

        if (index<N-2){
        index++;
        for(i=index-1, aux=0; i<index+2; i++,aux++){
          for (j=0;j<M;j++){
            enviada[(aux*M)+j]=w[i][j];
          }
        }
        //printf("Enviei a linha %d  da Iteracao %d \n", index, iter );       
        acabado=0;
        MPI_Send(&acabado, 1, MPI_INT,status.MPI_SOURCE,index,MPI_COMM_WORLD);
        MPI_Send(enviada, 3*M, MPI_DOUBLE,status.MPI_SOURCE,index,MPI_COMM_WORLD);

      }
      }



      for (env=1; env <numtasks; env++){
        acabado=1;
        MPI_Send(&acabado, 1, MPI_INT,env,env,MPI_COMM_WORLD);
      }



    for(i=0;i<N;i++)
      for(j=0;j<M;j++){
                temp=w[i][j];
                w[i][j]=u[i][j];
                u[i][j]=temp;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux ? diff : aux;
                printf(" diff: %lf \n",diff);
         }

        printf(" diff 2 : %lf \n",diff);

    }
    for (env=1; env <numtasks; env++){
        acabado=2;
        MPI_Send(&acabado, 1, MPI_INT,env,env,MPI_COMM_WORLD);
        printf(" HERE %d \n",iter);
      }
            

   }  
 
   /* End of master code */

   /************************* workers code **********************************/
  else{
      /* Receive rows from master */
      // Recebe indice
      //MPI_Recv(&index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); 
    while(1){
        while(1){
            MPI_Recv(&acabado, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
            if (acabado==0){
            MPI_Recv(enviada, 3*M, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status); 
            //
                //for(j=M;j<2*M;j++){
                //    printf("ANTES DO COMPUTE: linha/enviada[%d]=%lf\n",j,enviada[j]);
                //}
                //printf("antes do compute\n");
                compute(enviada, linha,3);
                //for(j=0;j<N;j++){
                    //printf("DEPOIS DO COMPUTE:linha[%d]=%lf\n",j,linha[j]);
                //}
                //
                //Envia computação para o Master
                //MPI_Send(&index,1, MPI_INT,0,0,MPI_COMM_WORLD); //check
                MPI_Send(linha,M, MPI_DOUBLE,0,status.MPI_TAG,MPI_COMM_WORLD); //check
                //printf("B\n");
            }
            else if(acabado==2){
                printf("Entrei no salto \n");
                goto end;       
            }
            else{
                printf("Entrei no break \n");
                break;
            }
    }
  }
}
    end:
      MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
      time =  MPI_Wtime() - time;
      printf("Elapsed time %lf\n",time);
      MPI_Finalize();
}


