#define _GNU_SOURCE
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/param.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define NUMSTEPS 1000000

int main() {
    struct timespec start, end;
    int comm_sz, my_rank;
    double pi = 0.0, step, local_sum = 0.0;
    long long i;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
    step = 1.0 / (double) NUMSTEPS;

    int workers = comm_sz - 1;

    clock_gettime(CLOCK_MONOTONIC, &start);

    if (my_rank != 0) {
      
        int base = NUMSTEPS / workers;
        int remainder = NUMSTEPS % workers;

        int my_count = base + (my_rank <= remainder ? 1 : 0);
        
        int my_start;
        if (my_rank <= remainder)
            my_start = (my_rank - 1) * (base + 1);
        else
            my_start = remainder * (base + 1) + (my_rank - remainder - 1) * base;

        int my_end = my_start + my_count;
       
        for (i = my_start; i < my_end; i++) {
            double x = (i + 0.5) * step;
            local_sum += 4.0 / (1.0 + x * x);
        }

        MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    } else {
        
        double recv_val;
        for (int src = 1; src < comm_sz; src++) {
            MPI_Recv(&recv_val, 1, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            pi += recv_val;
        }

        pi *= step;

        clock_gettime(CLOCK_MONOTONIC, &end);
        u_int64_t diff = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;

        printf("PI is %.20f\n", pi);
        printf("elapsed time = %llu nanoseconds\n", (long long unsigned int) diff);
    }

    MPI_Finalize();
    return 0;
}
