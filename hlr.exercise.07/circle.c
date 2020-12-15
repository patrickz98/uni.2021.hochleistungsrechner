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
    srand(28051998 * process);
    
    int bufferSize = (num / processes);
    int overflow = num % processes;

    // Check if number can be split evenly across all processes
    if (overflow != 0)
    {
        bufferSize++;
    }

    int generate = (num / processes);
    if (process < overflow)
    {
        // printf("[ %d ] --> Overflow\n", process);
        generate++;
    }

    printf("[ %d ] --> generate: %d\n", process, generate);

    int numbers[bufferSize];

    for (int inx = 0; inx < bufferSize; inx++)
    {
        numbers[ inx ] = (inx < generate)
            ? rand() % 25
            : 0;

        printf("[ %d ] --> numbers[ %d ] = %d\n", process, inx, numbers[ inx ]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("[ %d ] --> Process exit!\n", process);
    MPI_Finalize();
}