#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
extern void draw_heat(int nx, int ny);       /* X routine to create graph */

#define M      200                       /* x dimension of problem grid /Collumns*/
#define N      200                       /* y dimension of problem grid /Rows */
#define ITER       3               /* number of time steps */
#define MASTER      0                  /* taskid of first process */
#define BLOCOTAM  42
double **G1;
double **G2;



void fillMatrix(){
    int i,j;
        //CIMA
            for ( i = 0; i < M ; i++ ){
                G2[0][i] = 100;
                G1[0][i] = 100;
            }
        //BAIXO
            for ( i = 0; i < M ; i++ ){
                G2[N-1][i] = 100;
                G1[N-1][i] = 100;
            }
        //ESQUERDA
            for ( i = 1; i < N-1 ; i++ ){
                G2[i][0] = 100;
                G1[i][0] = 100;
            }
        //DIREITA    
            for ( i = 1; i < N -1 ; i++ ){
                G2[i][M-1] = 100;
                G1[i][M-1] = 100;
            }
        //RESTA    
            for ( i=1; i < N - 1 ; i++)
                for( j=1; j< M - 1 ; j++){
                    G2[i][j] = 0;
                }
    
}



/*****************************************************************************
 *  subroutine init
 *****************************************************************************/

void init(){
    int i;
    G1 = (double **) malloc(N*sizeof(double));
    G2 = (double **) malloc(N*sizeof(double));
    
    for(i = 0; i < N; i++){
        G1[i] = (double *) malloc(M*sizeof(double));
        G2[i] = (double *) malloc(M*sizeof(double));
    }
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
                 root at 0    root at 1
Left child        ix*2 + 1     ix*2
Right child       ix*2 + 2     ix*2 + 1
Parent            (ix-1)/2     ix/2
*/


int main (int argc, char *argv[]){

int taskid,                     /* this task's unique id */
  numtasks,                   /* number of tasks */
  source,               /* to - from for message send-receive */
  index,              /* loop variables */
  proc,
  start;
int iter;
MPI_Status status;
int j,i;
double ** temp;
int blocotam = BLOCOTAM;
double time;
double linha[M*(blocotam-2)];//, enviada pelos slaves;

double enviada [N*M];
int aux, nrslave;

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
      fillMatrix();


	 for(iter=0; iter<ITER;iter++){
        //printf("====Iteração nr %d===\n",iter);

        for(i=0, aux=0; i<N; i++,aux++){ //CONFIRMAR ESTE IF 
          for (j=0;j<M;j++){
            enviada[(aux*M)+j]=G1[i][j];
            //printf("Valor enviado G1[%d][%d]=%f com index %d \n", i,j,G1[i][j],index);
          }
        }


    MPI_Bcast(enviada, M*N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);




    int x = N/(blocotam-2);
//printf(" O VALOR DE X E %d\n", x);
for (; x>0; x--){
        
        MPI_Recv(linha, (M*(blocotam-2)), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        //Atualization of G2
        int line;
        line =(status.MPI_TAG-1)*(blocotam-2)+1;

        //printf(" ESTOU AQUI VALOR DE LINHA %d \n",line);
        for (i=0;i<(blocotam-2)*M ;i++){
         // printf(" X a %d valor de linha[%lf]\n",x,linha[i]);
        }
        for(i=0;i<(blocotam-2)*M && line < N-1;line++){
          for(j=0; j < M;j++,i++){
            G2[line][j]= linha[i];
            //printf("valor de G2[%d][%d] = %lf \n",line,j,G2[line][j]);
            //printf("Iteracao nr: %d G2[%d][%d]=%lf\n",iter,line,j,G2[line][j] );
        }
      }
  }
      

          for(i=0;i<N;i++)
      for(j=0;j<M;j++){
        G1[i][j]=G2[i][j];
        //  printf("G2[%d][%d]=%f \n",i,j,G1[i][j]);
         }
    }

}

   /************************* workers code **********************************/
  


  else{

    for (i=0;i<ITER;i++){


      MPI_Bcast(enviada, M*N, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

      	double computeover[blocotam*M];

      	/*for(i=0;i<N*M;i++){
      			printf("Envida %d = %f \n", i, enviada[i]);
      	}*/		
      	//if (taskid == 1){
      	start = (taskid-1)*((blocotam-2)*M);
      	/*}
        else{
      		start = (taskid-1)*(blocotam-2)*(M-1);
      	}*/

      	for (j=start, aux =0; j< start +((blocotam)*M) && start <N*M; j++,aux++){
      		computeover[aux]= enviada[j];
      		//printf("%d Compute %d = %f \n", taskid, aux, computeover[aux]);
      	}

        
        //for(j=0;j<N*(blocotam-1);j++){
            //printf("Nodo : %d ANTES DO COMPUTE:linha[%d]=%lf\n", taskid,j,computeover[j]);
        //}
        compute(computeover, linha, blocotam);


         //for(j=0;j<N*(blocotam-2);j++){
         //   printf("Nodo : %d DEPOIS DO COMPUTE:linha[%d]=%lf\n", taskid,j,linha[j]);
        //}

        MPI_Send(linha,M*(blocotam-2), MPI_DOUBLE,0,taskid,MPI_COMM_WORLD); //check
        //printf("B\n");
    
    }
  }
        MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
        time =  MPI_Wtime() - time;
        printf("Elapsed time %lf\n",time);
        MPI_Finalize();
}







