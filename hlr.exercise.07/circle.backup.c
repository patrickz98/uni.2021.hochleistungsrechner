#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>

int main(int argc, char* argv[])
{
    //
    // Arg parsing shit
    //

    if (argc != 2)
    {
        fprintf(stderr, "Only one parameter!\n");
        return -1;
    }

    int num = atoi(argv[1]);
    printf("N=%d\n", num);

    //
    // MPI setup
    //

    MPI_Init(NULL, NULL);

    int processes;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    int process;
    MPI_Comm_rank(MPI_COMM_WORLD, &process);

    // Init random with seed
    // srand(time(NULL));
    srand(28051998);
    
    int sliceChunks = (num / processes);
    int sliceStart  = (process       * sliceChunks);
    int sliceEnd    = ((process + 1) * sliceChunks);

    // Fix rounding errors
    if ((process + 1) == processes)
    {
        sliceEnd = num;
    }

    int sliceSize = sliceEnd - sliceStart;

    printf("[ %d ] --> sliceChunks=%d sliceStart=%d, sliceEnd=%d\n",
        process, sliceChunks, sliceStart, sliceEnd);

    int *numbers = malloc(sliceSize * sizeof(int));

    for (int inx = 0; inx < sliceSize; inx++)
    {
        numbers[ inx ] = rand() % 25;
        printf("[ %d ] --> numbers[ %d ] = %d\n", process, (sliceStart + inx), numbers[ inx ]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("[ %d ] --> Process exit!\n", process);
    MPI_Finalize();
}