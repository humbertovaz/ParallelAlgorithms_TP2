#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <omp.h>


#define N 200
#define TOL tol
#define CLUSTER 100
#define TIME_RESOLUTION 1000000 /* time measuring resolution (us) */
#define epsilon FLT_EPSILON

#define NB 8 /*block size */
#define min(a,b) ( ((a) < (b)) ? (a) : (b) )






long long unsigned initial_time;
struct timeval t;

double u[N][N], w[N][N];
double tol = 1.0f / (N * N);


double uArray[N*N];
 double utmp [N*N];

void start(void) {
    gettimeofday(&t, NULL);
    initial_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;
}

long long unsigned stop(void) {
    gettimeofday(&t, NULL);
    long long unsigned final_time = t.tv_sec * TIME_RESOLUTION + t.tv_usec;

    return final_time - initial_time;
}




/*************** Inicializaçao *****************
% SET BOUNDARY VALUES 
% Temperature is zero at top and 100 on the other boundaries
u(1,1:N)=100;  % lower boundary 
u(1:N,1)=100;  % left boundary 
u(1:N,N)=100;  % right boundary 
u(N,1:N)=0;    % boundary above

w(1,1:N)=100;  % lower boundary 
w(1:N,1)=100;  % left boundary 
w(1:N,N)=100;  % right boundary 
w(N,1:N)=0;    % boundary above

% initial values for interior points
u(2:N-1,2:N-1)=50;

***********************************************/



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

void updateMatrix(){
      memcpy(u,w, N*N*sizeof(double));
        /*
        Talvez usar swap de apontadores
        temp = u;
        w = u;
        u = temp;
        */

      
}



/************************************** Parallel ****************************************/


int poissonGS_parallel(){
    double diff = TOL +1, aux, diffmax;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
        #pragma omp parallel for reduction(max : diff) shared(w, u) private(i,j, aux)
        for (i = 1; i < (N-1); i++)
            for (j = (i %2 == 0? 2 : 1); j < (N-1); j+=2){
                w[i][j] = ( u[i][j+1] + u[i+1][j] +u[i-1][j] + u[i][j-1] ) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux?  : aux;
                //printf ("diff = %lf \n",diff);
            }
        diffmax = diff;

        #pragma omp parallel for reduction(max : diff) shared(w, u) private(i,j, aux)
        for (i = 1; i < (N-1); i++)
            for (j = (i %2 == 0? 1 : 2); j < (N-1); j+=2){
                u[i][j] = w[i][j];
                w[i][j] = ( w[i][j+1] + w[i+1][j] +w[i-1][j] + w[i][j-1] ) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux?  : aux;
                //printf ("diff2 = %lf \n",diff);
            }    
        diff = diff > diffmax ? diff : diffmax;
        updateMatrix();
        n_iter++;           

    }
    return n_iter;
}






int poissonNaive_parallel(){
    double diff = TOL +1, aux;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
        #pragma omp parallel for reduction(max : diff) shared(w, u) private(i,j, aux)
        for (i = 1; i < (N-1); i++)
            for (j = 1; j < (N-1); j++){
                w[i][j] = (u[i-1][j] + u[i][j-1] + u[i][j+1] + u[i+1][j]) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux? diff : aux;
            }
        updateMatrix();
        n_iter++;           

    }
    return n_iter;
}






int poissonGS_SOR2_parrallel(){
    double diff = TOL +1, aux, updateval;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
        #pragma omp parallel for reduction(max : diff) shared(w, u) private(i,j, aux)
        for (i = 1; i < (N-1); i++)
            for (j = 1; j < (N-1); j++){
                updateval = (( u[i][j+1] + u[i+1][j] +w[i-1][j] + w[i][j-1] ) *0.25) - w[i][j];
                w[i][j] += updateval * (1+sin (M_PI/(N+1)));
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux? diff : aux;
            }
        updateMatrix();
        n_iter++;           

    }
    return n_iter;
}







/*********************************** Parallel blocks ************************************/


/*
 * Blocked Jacobi solver: one iteration step
 */

/*
double relax_jacobi ()
{
	double diff, aux=0.0;
	int nbx, bx, nby, by,ii,jj,i,j;
    int n_iter=0;
	nbx = omp_get_max_threads();//NB;
	bx = N/nbx + ((N%nbx) ? 1 : 0);//sizex/nbx;
	nby = 1;//NB;
	by = N/nby;
    diff = TOL +1;
    while (diff > TOL){
	#pragma omp parallel for reduction(max : diff) shared(w, u) private(i,ii,j,jj, aux)
	for (ii=0; ii<nbx; ii++) {
		for ( jj=0; jj<nby; jj++)  {
			for (i=1+ii*bx; i<=min((ii+1)*bx, N-2); i++) {
				for (j=1+jj*by; j<=min((jj+1)*by, N-2); j++) {
                    w[i][j] = ( u[i][j+1] + u[i+1][j] +u[i-1][j] + u[i][j-1] ) *0.25;
					aux = (w[i][j] - u[i][j]);
                    diff = diff > aux? diff : aux; 
                    printf ("diff2 = %lf \n",diff);
				}
			}
		}
	}

        updateMatrix();
        n_iter++;           
    }
    return n_iter;
}






double relax_jacobi ()
{
    int ii, jj, i, j;
    int n_iter=0;
    double diff = TOL +1, aux;
    int nbx, bx, nby, by;
    nbx = NB;
    bx = N/nbx;
    nby = NB;
    by = N/nby;
    while (diff > TOL){
        diff = 0;
        #pragma omp parallel for reduction(max : diff) shared(w, u) private(i,ii,j,jj, aux)
        for (ii=0; ii<nbx; ii++)
            for (jj=0; jj<nby; jj++) 
                for (i=1+ii*bx; i<=min((ii+1)*bx, N-2); i++) 
                    for (j=1+jj*by; j<=min((jj+1)*by, N-2); j++) {
                        w[i][j] = ( u[i][j+1] + u[i+1][j] +u[i-1][j] + u[i][j-1] ) *0.25;
                        aux = (w[i][j] - u[i][j]);
                        diff = diff > aux? diff : aux;
                        printf ("diff2 = %lf \n",diff);
                }
        updateMatrix();
        n_iter++;           
    }
    return n_iter;
}


*/



/******************* DRAW *********************/


void printi(char* name){
	double value;
	FILE *f = fopen(strcat(name,".ppm"), "wb");
	fprintf(f, "P6\n%i %i 255\n", N, N);
	for (int y=0; y<N; y++) {
	    for (int x=0; x<N; x++) {
	    	value=w[y][x];
	    	
	    	if(value<50){
		        fputc(0, f);   // 0 .. 255
		        fputc(value*(255/50), f);   // 0 .. 255
	        	fputc(255-(value*(255/50)), f);   // 0 .. 255
	    	}else{
	    		fputc((value-50)*(255/50), f);   // 0 .. 255
		        fputc(255-((value-50)*(255/50)), f);   // 0 .. 255
	        	fputc(0, f);   // 0 .. 255
	    	}
	    }
	}
	fclose(f);
}





int main(int argc, char** argv){
    int n_iter;
    long long unsigned times[3];
    init();
    start();
    n_iter= poissonNaive_parallel(); 
    times[0] = stop();
    printf("Parallel Jacobi. Number of iterations: %d. Time elapsed: %llu usecs\n", n_iter, times[0]);
   
    char outputname[100];
    strcpy(outputname,"JacobiPAR");
	printi(outputname);



    init();
    start();
    n_iter= poissonGS_parallel();
    times[1] = stop();
    printf("Parallel poissonGS (RedBlack). Number of iterations: %d. Time elapsed: %llu usecs\n", n_iter, times[1]);


    strcpy(outputname,"PGSPAR");
	printi(outputname);





    init();
    start();
    n_iter= poissonGS_SOR2_parrallel();
    times[2] = stop();
    printf("Parallel PoissonGS_SOR2. Number of iterations: %d. Time elapsed: %llu usecs\n", n_iter, times[2]);


    strcpy(outputname,"PGSSOR2PAR");
	printi(outputname);




    /*Blocked version
    init();
    start();
    n_iter= relax_jacobi();
    times[3] = stop();
    printf("Parallel Jacobi Blocks. Number of iterations: %d. Time elapsed: %llu usecs\n", n_iter, times[3]);

    strcpy(outputname,"JacobiBlocksPAR");
	printi(outputname);

    */
}