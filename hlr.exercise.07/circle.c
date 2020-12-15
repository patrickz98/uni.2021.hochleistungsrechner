#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>

#define TAG_START_NUMBER 0
#define TAG_SEND_NUMBERS 1
#define TAG_EXIT         2
#define TAG_NUMBER_COUNT 3

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
    // printf("N=%d\n", num);

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
    int overflow   = (num % processes);

    // Check if number can be split evenly across all processes
    if (overflow != 0)
    {
        bufferSize++;
    }

    int numberCount = (num / processes);
    if (process < overflow)
    {
        numberCount++;
    }

    int numbers[bufferSize];

    for (int inx = 0; inx < bufferSize; inx++)
    {
        numbers[ inx ] = (inx < numberCount)
            ? rand() % 25
            : 0;
    }

    for (int inx = 0; inx < numberCount; inx++)
    {
        printf("[ %d ] --> numbers[ %d ] = %d\n", process, inx, numbers[ inx ]);
    }

    //
    // MPI Stuff
    //

    MPI_Status status;

    if (process == 0)
    {
        //
        // Send numbers[0] to last process
        //

        int destination = (processes - 1);
        MPI_Send((void*) &numbers[0], 1, MPI_INT, destination, TAG_START_NUMBER, MPI_COMM_WORLD);
        printf("[ %d ] --> Send: %d\n", process, numbers[0]);
    }

    int stopNumber = -1;
    int brodcastRoot = (processes - 1);
    int isProc0 = (process == brodcastRoot);

    if (isProc0)
    {
        //
        // Receive stop nummber from process[0]
        //

        MPI_Recv(&stopNumber, 1, MPI_INT, 0, TAG_START_NUMBER, MPI_COMM_WORLD, &status);
        printf("[ %d ] --> Recv: %d\n", process, stopNumber);
    }

    // Process source nummber
    int source = (process + processes - 1) % processes;

    // Process destination nummber
    int destination = (process + processes + 1) % processes;

    // printf("[ %d ] --> source: %d, destination: %d\n", process, source, destination);

    int exit = 0;
    
    while (exit == 0)
    {
        MPI_Send((void*) &numbers, bufferSize, MPI_INT, destination, TAG_SEND_NUMBERS, MPI_COMM_WORLD);
        MPI_Send((void*) &numberCount, 1, MPI_INT, destination, TAG_NUMBER_COUNT, MPI_COMM_WORLD);

        MPI_Recv(&numbers, bufferSize, MPI_INT, source, TAG_SEND_NUMBERS, MPI_COMM_WORLD, &status);
        MPI_Recv(&numberCount, 1, MPI_INT, source, TAG_NUMBER_COUNT, MPI_COMM_WORLD, &status);

        if (isProc0)
        {
            exit = (stopNumber == numbers[0]);
        }

        MPI_Bcast((void*) &exit, 1, MPI_INT, brodcastRoot, MPI_COMM_WORLD);
    }

    for (int inx = 0; inx < numberCount; inx++)
    {
        printf("[ %d ] --> final: numbers[ %d ] = %d\n", process, inx, numbers[ inx ]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("[ %d ] --> Process exit!\n", process);
    MPI_Finalize();
}