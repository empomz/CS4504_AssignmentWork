#include <stdio.h>
#include <omp.h>
#include <time.h>
#include <stdlib.h>

#define N 2048
#define FactorIntToDouble 1.1;

double firstMatrix [N] [N] = {0.0};
double secondMatrix [N] [N] = {0.0};
double matrixMultiResult [N] [N] = {0.0};


void matrixMulti()
{
    for(int row = 0 ; row < N ; row++){
        for(int col = 0; col < N ; col++){
            double resultValue = 0;
            for(int transNumber = 0 ; transNumber < N ; transNumber++) {
                resultValue += firstMatrix [row] [transNumber] * secondMatrix [transNumber] [col] ;
            }

            matrixMultiResult [row] [col] = resultValue;
        }
    }
}


void matrixInit()
{
    for(int row = 0 ; row < N ; row++ ) {
        for(int col = 0 ; col < N ;col++){
            srand(row+col);
            firstMatrix [row] [col] = ( rand() % 10 ) * FactorIntToDouble;
            secondMatrix [row] [col] = ( rand() % 10 ) * FactorIntToDouble;
        }
    }
}



int main(int argc, char *argv[])
{
    int thread_count = strtol(argv[1], NULL, 10);

#pragma omp parallel num_threads(thread_count)
    matrixInit();

    //clock_t t1 = clock();
    //matrixMulti();
    //clock_t t2 = clock();
    //printf("time: %ld", t2-t1);

    double t1 = omp_get_wtime();
    matrixMulti();
    double t2 = omp_get_wtime();
    printf("serial   time: %3f\n", ((double)t2 - t1) / CLOCKS_PER_SEC * 1000000.0);

    return 0;
}
