#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <mpi.h>

int myrank; // Current rank
int ranksize; // Number of processes in communicator

int main(int argc, char **argv)
{
    int error;
    error = MPI_Init(&argc, &argv);
    if (error != MPI_SUCCESS) {
        fprintf(stderr, "ERROR: CAN'T MPI INIT\n");
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
    MPI_Barrier(MPI_COMM_WORLD);



    printf("Hello\n");
    MPI_Finalize();
    return 0;
}