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
#define TAG_PRINT_SYNC   4

void printSynced(int proc, int processes, int *numbers, int bufferSize)
{
    MPI_Status status;

    if (proc == 0)
    {
        for (int inx = 0; inx < bufferSize; inx++)
        {
            int num = numbers[inx];

            if (num >= 0)
            {
                printf("%2d ", num);
            }
            else
            {
                printf("   ");
            }
        }

        printf(" | ");

        int l = 1;
        MPI_Send((void*) &l, 1, MPI_INT, 1, TAG_PRINT_SYNC, MPI_COMM_WORLD);

        return;
    }

    int l = 0;
    MPI_Recv((void*) &l, 1, MPI_INT, (proc - 1), TAG_PRINT_SYNC, MPI_COMM_WORLD, &status);

    for (int inx = 0; inx < bufferSize; inx++)
    {
        int num = numbers[inx];

        if (num >= 0)
        {
            printf("%2d ", num);
        }
        else
        {
            printf("   ");
        }
    }

    if (proc != (processes - 1))
    {
        printf(" | ");
        MPI_Send((void*) &l, 1, MPI_INT, (proc + 1), TAG_PRINT_SYNC, MPI_COMM_WORLD);
    }
    else
    {
        printf("\n");
    }
}

void circle(int num)
{
    //
    // MPI setup
    //

    MPI_Init(NULL, NULL);

    MPI_Status status;

    int processes;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);

    int process;
    MPI_Comm_rank(MPI_COMM_WORLD, &process);
    
    //
    // Set a static buffer size across all process
    //

    int bufferSize = (num / processes);
    int overflow   = (num % processes);

    // Check if number can be split evenly across all processes
    // if not increase buffer size
    if (overflow != 0)
    {
        bufferSize++;
    }

    // Calculate how many random numbers need to be generated
    int numberCount = (num / processes);
    if (process < overflow)
    {
        numberCount++;
    }

    // Put numbers on the heap
    int *numbers = malloc(sizeof(int) * bufferSize);

    // Init "random" with seed
    srand(28051998 * process);

    // Fill numbers array with random numbers
    // The -1 signals overflow
    for (int inx = 0; inx < bufferSize; inx++)
    {
        numbers[ inx ] = (inx < numberCount)
            ? rand() % 25
            : -1;
    }

    // Check for only one process
    if (processes <= 1)
    {
        printf("Array befor and after rotate: ");

        for (int inx = 0; inx < bufferSize; inx++)
        {
            int number = numbers[ inx ];

            if (number >= 0)
            {
                printf("%2d ", number);
            }
        }

        printf("\n");

        MPI_Finalize();
        free(numbers);

        return;
    }

    // Print out initial array
    if (process == 0) printf("Array befor rotate: ");
    printSynced(process, processes, numbers, bufferSize);

    //
    // Send stop nummber from P0 to last process
    //

    if (process == 0)
    {
        //
        // Send numbers[0] to last process
        //

        int destination = (processes - 1);
        MPI_Send((void*) &numbers[0], 1, MPI_INT, destination, TAG_START_NUMBER, MPI_COMM_WORLD);
        // printf("[ %d ] --> Send: %d\n", process, numbers[0]);
    }

    int stopNumber = -1;
    int brodcastRoot = (processes - 1);
    int isLastProc = (process == brodcastRoot);

    if (isLastProc)
    {
        //
        // Receive stop nummber from process[0]
        //

        MPI_Recv(&stopNumber, 1, MPI_INT, 0, TAG_START_NUMBER, MPI_COMM_WORLD, &status);
        // printf("[ %d ] --> Recv: %d\n", process, stopNumber);
    }

    // Calculate source and destination process number for message exchange
    int source      = (process + processes - 1) % processes;
    int destination = (process + processes + 1) % processes;

    //
    // Send and receive number array
    //

    int iterations = 0;
    int exit = 0;
    while (exit == 0)
    {
        MPI_Send((void*) numbers, bufferSize, MPI_INT, destination, TAG_SEND_NUMBERS, MPI_COMM_WORLD);
        MPI_Recv(        numbers, bufferSize, MPI_INT,      source, TAG_SEND_NUMBERS, MPI_COMM_WORLD, &status);

        if (isLastProc)
        {
            // Check if number matches number from first process
            exit = (stopNumber == numbers[0]);
        }

        // Brodcast termination
        MPI_Bcast((void*) &exit, 1, MPI_INT, brodcastRoot, MPI_COMM_WORLD);
        iterations++;
    }

    //
    // Print final array after exchange
    //

    if (process == 0) printf("Array after rotate: ");
    printSynced(process, processes, numbers, bufferSize);

    //
    // MPI clean-up
    //

    MPI_Barrier(MPI_COMM_WORLD);

    if (process == 0) printf("Iterations: %d\n", iterations);
    free(numbers);

    MPI_Finalize();
}

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

    circle(num);
}