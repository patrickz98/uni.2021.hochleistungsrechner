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

    // Threads 1 bis n sollen bestimmte Anweisungen ausführen
    if (world_rank != 0)
    {
        //
        // 1. Wir senden den Hostnamen.
        //

        MPI_Send((void*) processor_name, name_len, MPI_CHAR, 0, 1, MPI_COMM_WORLD);

        //
        // 2. Wir senden die lokale Zeit als String.
        //

        // Wir speichern die lokale Zeit als String
        // im Format Jahr-Monat-Tag Stunde-Minute-Sekunde
        time_t now = time(NULL);
        struct tm tm_now;

        // Wir messen die Zeit bis auf die Sekunde genau
        localtime_r(&now, &tm_now);

        // Wir speichern die Zeit als String
        char time_string[100];
        strftime(time_string, sizeof(time_string), "%Y-%m-%d %H:%M:%S", &tm_now);

        MPI_Send((void*) time_string, 100, MPI_CHAR, 0, 2, MPI_COMM_WORLD);

        // 
        // 3. Wir senden die Mikrosekunden
        //

        // Zum Abspeichern der Mikrosekunden seit der letzten vollen Sekunde
        struct timeval current_time;
        gettimeofday(&current_time, NULL);
        int microsec = current_time.tv_usec;

        MPI_Send((void*) &microsec, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
    }
    else
    {
        //
        // Anweisungen für Thread 0
        //

        // Wir definieren Arrays um die verschiedenen Hostnamen und Timestamps zu speichern
        char* hostname[world_size];
        char* time_string[world_size];
        int microsec[world_size];

        for (int inx = 1; inx < world_size; inx++)
        {
            //
            // 1. Wir lesen den Hostnamen aus und speichern ihn im Array
            //

            // Wir nutzen MPI_Probe um die Größe der Nachricht abzuschätzen
            MPI_Probe(inx, 1, MPI_COMM_WORLD, &status);

            int length_name;
            MPI_Get_count(&status, MPI_CHAR, &length_name);
            char* buf_name = malloc(sizeof(char) * length_name);
            MPI_Recv(buf_name, length_name, MPI_CHAR, inx, 1, MPI_COMM_WORLD, &status);
            hostname[inx] = buf_name;

            //
            // 2. Wir lesen den Timestamp aus und speichern ihn im Array
            //

            // Wir nutzen MPI_Probe um die Größe der Nachricht abzuschätzen
            MPI_Probe(inx, 2, MPI_COMM_WORLD, &status);
            
            int length_time;
            MPI_Get_count(&status, MPI_CHAR, &length_time);
            char* buf_time = malloc(sizeof(char) * length_time);

            MPI_Recv(buf_time, length_time, MPI_CHAR, inx, 2, MPI_COMM_WORLD, &status);
            time_string[inx] = buf_time;

            //
            // 3. Wir lesen die Mikrosekunden aus und speichern sie im Array
            //
            
            MPI_Recv(&microsec[inx], 1, MPI_INT, inx, 3, MPI_COMM_WORLD, &status);
        }

        // Für einen einzelnen Thread brauchen wir keine Ausgabe
        if (world_size > 1)
        {
            for (int inx = 1; inx < world_size; inx++)
            {
                // Thread 0 fängt mit der Ausgabe an
                printf("%s: %s.%d\n", hostname[inx], time_string[inx], microsec[inx]);
            }
        }
    }

    // Jeder Threads warten hier, bis alle Threads an dieser Stelle angekommen sind. Erst dann werden die Threads beendet.
    MPI_Barrier(MPI_COMM_WORLD);

    // DIREKT vor dem Beenden kommt die exit-Ausgabe. Aufgabenstellung war etwas unklar,
    // ob die Barriere vielleicht nach dieser Ausgabe erfolgen sollte.
    printf("Rang %d beendet jetzt!\n", world_rank);

    // Finalize the MPI environment.
    MPI_Finalize();
}
