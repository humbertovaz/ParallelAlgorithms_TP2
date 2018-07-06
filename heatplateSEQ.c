#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>


#define N 2000
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

      /*Talvez usar swap de apontadores*/
}

/*
*Tradução ex do dado
*/
int Jacobi(){
    double diff = TOL +1, aux;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
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



int poissonGS(){
    double diff = TOL +1, aux;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
        for (i = 1; i < (N-1); i++)
            for (j = 1; j < (N-1); j++){
                //w[i][j] = (w[i-1][j] + w[i][j-1] + u[i][j+1] + u[i+1][j]) *0.25;
                w[i][j] = ( u[i][j+1] + u[i+1][j] +w[i-1][j] + w[i][j-1] ) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux? diff : aux;
            }
        updateMatrix();
        n_iter++;           

    }
    return n_iter;
}



int poissonGS2(){
    double diff = TOL +1, aux,diffmax;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;

        for (i = 1; i < (N-1); i++)
            for (j = (i %2 == 0? 2 : 1); j < (N-1); j+=2){
                w[i][j] = ( u[i][j+1] + u[i+1][j] +u[i-1][j] + u[i][j-1] ) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux ?  : aux;
            }
        diffmax = diff;
        for (i = 1; i < (N-1); i++)
            for (j = (i %2 == 0? 1 : 2); j < (N-1); j+=2){
                u[i][j]=w[i][j];
                w[i][j] = ( w[i][j+1] + w[i+1][j] +w[i-1][j] + w[i][j-1] ) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux? : aux;
            }    
        diff = diff > diffmax ? diff : diffmax;
        updateMatrix();
        n_iter++;    

        for (i = 1; i < (N-1); i++)
            for (j = 1; j < (N-1); j++){
                
                w[i][j] = ( u[i][j+1] + u[i+1][j] +w[i-1][j] + w[i][j-1] ) *0.25;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux? diff : aux;
            }
        updateMatrix();
        n_iter++;           

    }
    return n_iter;
}



/*
* successive overrelaxation (SOR) com valor estatico 1.5
*/


int poissonGS_SOR(){
    double diff = TOL +1, aux, updateval;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
        for (i = 1; i < (N-1); i++)
            for (j = 1; j < (N-1); j++){
                //w[i][j] = (w[i-1][j] + w[i][j-1] + u[i][j+1] + u[i+1][j]) / 4.0;
                updateval = (( u[i][j+1] + u[i+1][j] +w[i-1][j] + w[i][j-1] ) *0.25) - w[i][j];
                w[i][j] += updateval * 1.5;
                aux = (w[i][j] - u[i][j]);
                diff = diff > aux? diff : aux;
            }
        updateMatrix();
        n_iter++;           

    }
    return n_iter;
}




/*
* successive overrelaxation (SOR) com valor estatico adaptado (1+sin (M_PI/(N+1)))<2
*/


int poissonGS_SOR2(){
    double diff = TOL +1, aux, updateval;
    int n_iter = 0;
    int i, j;
    while (diff > TOL){
        diff = 0;
        for (i = 1; i < (N-1); i++)
            for (j = 1; j < (N-1); j++){
                //w[i][j] = (w[i-1][j] + w[i][j-1] + u[i][j+1] + u[i+1][j]) / 4.0;
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







/************** SEQ with blocks ********************/

/*
 * Blocked Jacobi solver:
 */
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
        for (ii=0; ii<nbx; ii++)
            for (jj=0; jj<nby; jj++) 
                for (i=1+ii*bx; i<=min((ii+1)*bx, N-2); i++) 
                    for (j=1+jj*by; j<=min((jj+1)*by, N-2); j++) {
                        w[i][j] = ( u[i][j+1] + u[i+1][j] +u[i-1][j] + u[i][j-1] ) *0.25;
                        aux = (w[i][j] - u[i][j]);
                        diff = diff > aux? diff : aux;
                }
        updateMatrix();
        n_iter++;           
    }
    return n_iter;
}






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




/**************************************************************************************/



void debug_visualize_matrix(){
    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++) printf("%.1f ", u[i][j]);
        putchar('\n');
    }
}

int main(int argc, char** argv){
    int n_iter;
    long long unsigned times[6];

    init();
    start();
    n_iter= Jacobi();
    times[0] = stop();
    printf("Jacobi: Number of iterations: %d. Time elapsed: %llu usecs\n\n", n_iter, times[0]);

    char outputname[100];
    strcpy(outputname,"Jacobi");
	printi(outputname);


    init();
    start();
    n_iter= poissonGS();
    times[1] = stop();
    printf("PoissonGS: Number of iterations: %d. Time elapsed: %llu usecs\n\n", n_iter, times[1]);
    

    strcpy(outputname,"PGS1");
	printi(outputname);


    init();
    start();
    n_iter= poissonGS_SOR();
    times[2] = stop();
    printf("PoissonGS_SOR: Number of iterations: %d. Time elapsed: %llu usecs\n\n", n_iter, times[2]);
    


    strcpy(outputname,"PGSSOR");
	printi(outputname);



    init();
    start();
    n_iter= relax_jacobi();
    times[3] = stop();
    printf("Jacobi Blocks: Number of iterations: %d. Time elapsed: %llu usecs\n\n", n_iter, times[3]);
        

    strcpy(outputname,"JacobiBlocks");
	printi(outputname);

    
    init();
    start();
    n_iter= poissonGS2();
    times[4] = stop();
    printf("PoissonGS 2 : Number of iterations: %d. Time elapsed: %llu usecs\n\n", n_iter, times[4]);



    strcpy(outputname,"PGS2");
	printi(outputname);



    init();
    start();
    n_iter= poissonGS_SOR2();
    times[5] = stop();
    printf("PoissonGS_SOR 2 : Number of iterations: %d. Time elapsed: %llu usecs\n\n", n_iter, times[5]);


    strcpy(outputname,"PGSSOR2");
	printi(outputname);


    //debug_visualize_matrix();
    /*
    init();
    start();
    poissonGS_RB();
    printf("Solution found through poissonGS_RB. Time elapsed: %llu usecs\n\n", stop());

    //TODO: try different nthreads

    init();
    start();
    poissonGS_RB_SOR();
    printf("Solution found through poissonGS_RB_SOR. Time elapsed: %llu usecs\n\n", stop());
    */
    return 0;

}