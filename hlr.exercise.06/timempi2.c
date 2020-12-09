#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>

int main()
{
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor/host
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    MPI_Status status;

    if(world_rank != 0) /*Threads 1 bis n sollen bestimmte Anweisungen ausführen*/
    {
        /*Zum Abspeichern der Mikrosekunden seit der letzten vollen Sekunde*/
        struct timeval current_time;

        /*Wir speichern die lokale Zeit als String im Format Jahr-Monat-Tag Stunde-Minute-Sekunde*/
        time_t now = time(NULL) ;
        struct tm tm_now ;
        localtime_r(&now, &tm_now); //Wir messen die Zeit bis auf die Sekunde genau

        gettimeofday(&current_time, NULL);
        int microsec = current_time.tv_usec;//Wir messen die Zeit bis auf die Mikrosekunde genau

        char time_string[100];
        strftime(time_string, sizeof(time_string), "%Y-%m-%d %H:%M:%S", &tm_now);//Wir speichern die Zeit als String

        /*Wir senden den Hostnamen.*/
        MPI_Send((void*) processor_name, name_len, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
        /*Wir senden die lokale Zeit als String*/
        MPI_Send((void*) time_string, 100, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
        /*Wir senden die Mikrosekunden */
        MPI_Send((void*) &microsec, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
        //printf("sent from %d: %s\n", world_rank, processor_name);
    }
    else /*Anweisungen für Thread 0*/
    {
        /*Wir definieren Arrays um die verschiedenen Hostnamen und Timestamps zu speichern*/
        char* hostname[world_size];
        char* time_string[world_size];
        int microsec[world_size];

        for(int i = 1; i < world_size; i++)
        {
            /*Wir nutzen MPI_Probe um die Größe der Nachricht abzuschätzen*/
            MPI_Probe(i, 1, MPI_COMM_WORLD, &status);
            int length_name;
            MPI_Get_count(&status, MPI_CHAR, &length_name);
            char* buf_name = malloc(sizeof(char) * (length_name));
            /*Wir lesen den Hostnamen aus und speichern ihn im Array*/
            MPI_Recv(buf_name, length_name, MPI_CHAR, i, 1, MPI_COMM_WORLD, &status);
            hostname[i] = buf_name;

            /*Wir nutzen MPI_Probe um die Größe der Nachricht abzuschätzen*/
            MPI_Probe(i, 2, MPI_COMM_WORLD, &status);
            int length_time;
            MPI_Get_count(&status, MPI_CHAR, &length_time);
            char* buf_time = malloc(sizeof(char) * (length_time));
            /*Wir lesen den Timestamp aus und speichern ihn im Array*/
            MPI_Recv(buf_time, length_time, MPI_CHAR, i, 2, MPI_COMM_WORLD, &status);
            time_string[i] = buf_time;

            /*Wir lesen die Mikrosekunden aus und speichern sie im Array*/
            MPI_Recv(&microsec[i], 1, MPI_INT, i, 3, MPI_COMM_WORLD, &status);
        }

        if(world_size > 1) /*Für einen einzelnen Thread brauchen wir keine Ausgabe*/
        {
            for(int i = 1; i < world_size; i++)
            { /*Thread 0 fängt mit der Ausgabe an*/
                //printf("recieved from thread %d:\n", i);
                printf("%s : %s:%d\n", hostname[i], time_string[i], microsec[i]); //
            }

            /*Wir suchen die kleinste Anzahl an Mikrosekunden*/
            int min = INT_MAX;
            for(int i = 1; i < world_size; i++)
            {
                if(microsec[i] < min)
                {
                    min = microsec[i];
                }
            }
            printf("%d\n", min);
        }
    }

    /*Jeder Threads warten hier, bis alle Threads an dieser Stelle angekommen sind. Erst dann werden die Threads beendet.*/
    MPI_Barrier(MPI_COMM_WORLD);

    /*DIREKT vor dem Beenden kommt die exit-Ausgabe. Aufgabenstellung war etwas unklar, ob die Barriere vielleicht nach dieser Ausgabe erfolgen sollte.*/
    printf("Rang %d beendet jetzt!\n", world_rank);
    // Finalize the MPI environment.
    MPI_Finalize();
}
