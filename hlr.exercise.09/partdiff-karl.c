/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
	uint64_t  comm_size;      /* Anzahl der Prozesse die das Programm bearbeiten*/
	uint64_t  rang;           /* Rang des Prozesses                             */
	uint64_t  rows;           /* Anzahl der Zeilen die je ein Prozess bearbeitet*/
	int offset;               /* Differenz vom globalen zum lokalen Zeilenrang  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

static
int
calc_chunk (struct calculation_arguments* arguments)
{
    int error_flag = 0;

    uint64_t const N = arguments->N;

    int world_size = arguments->comm_size;

	/*Berechnung des chunks den je ein Prozess übernimmt. Entspricht hier der Anzahl an Zeilen der Matrix die vom jeweiligen Prozess berechnet werden.*/
	uint64_t rows = (N - 1) / world_size;

	/*Da sich die Matrix nur im Ausnahmefall glatt auf alle Prozesse aufteilen lässt, müssen einige Prozesse mehr Zeilen
	übernehmen als andere. In diesem Fall übernehmen die ersten x Prozesse eine Zeile mehr.*/
	if((N - 1) % world_size != 0)
    {
        if(arguments->rang < (N - 1) % world_size)
        {
            rows++;
        }
    }

    /*Jeder Prozess speichert 2 Randzeilen. Ab hier entspricht rows der vom Prozess gehaltenen Zeilen.*/
    rows += 2;

    /*offset Beschreibt die Position der jeweiligen Zeile in der Gesamtmatrix. Für die sin-Berechnung bei clalculate_parallel benötigt,
    sowie bei init_matrices_parallel und display_matrices_parallel.
    i + offset soll auf den globalen Rang der lokalen Zeile mit Index i zeigen*/
    int offset = 0;
    MPI_Status status;

	if(arguments->rang == 0)
    {
        arguments->offset = offset;

        offset += rows - 2; //Wir ziehen 2 Zeilen ab, eine für die Randzeile dieses Prozesses, eine für die Randzeile des nächsten Prozesses

        MPI_Send((void*) &offset, 1, MPI_INT, arguments->rang + 1, 1, MPI_COMM_WORLD);
    }
    else if(arguments->rang + 1 == arguments->comm_size)
    {
        MPI_Recv((void*) &offset, 1, MPI_INT, arguments->rang - 1, 1, MPI_COMM_WORLD, &status);

        arguments->offset = offset;
    }
    else
    {
        MPI_Recv((void*) &offset, 1, MPI_INT, arguments->rang - 1, 1, MPI_COMM_WORLD, &status);

        arguments->offset = offset;

        offset += rows - 2; //Wir ziehen 2 Zeilen ab, eine für die Randzeile dieses Prozesses, eine für die Randzeile des nächsten Prozesses

        MPI_Send((void*) &offset, 1, MPI_INT, arguments->rang + 1, 1, MPI_COMM_WORLD);
    }

    /*Debugging und testen*/
    //printf("rang: %" PRIu64 "\nrows: %" PRIu64 "\noffset: %d\n\n", arguments->rang, rows, arguments->offset);

    /*  Fehlerbehandlung
    Wenn ein Prozess nur die Randzeilen abspeichert, kommt es zu falschen Matrixausgaben. Wir brechen das Programm in diesem Fall ab.*/
    const int error = (rows <= 2) ? 1 : 0;

    /*Alle Prozesse erhalten den error-Wert. Wenn ein Prozess <= 2 Zeilen hält, ist die error_flag für alle Prozesse nun 1.*/
    MPI_Allreduce((void*) &error, (void*) &error_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    /*Prozess 0 gibt dem User Bescheid.*/
    if((arguments->rang == 0) && (error_flag != 0))
    {
        printf("Das Programm ist mit diesem Verhältnis von Interlines und Prozessen nicht ausführbar.\n");
        printf("Prozesse: %d\nInterlines: %" PRIu64 "\n", world_size, (N - 8) / 8);
        printf("Beende Programm...\n");
    }

    arguments->rows = rows;

    /*Falls ein Fehler aufgetaucht ist (bzw die Aufteilung zu Fehlern führt), gilt (error_flag != 0).
    Wir geben dies der main-Funktion zurück, welche das Programm beendet.*/
    return error_flag;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* allocateMatrices_parallel: allocates memory for a chunk of the matrices  */
/* ************************************************************************ */
static
void
allocateMatrices_parallel (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;
	uint64_t const rows = arguments->rows;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * rows * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
	    /*Speicherallokierung für die Zeilen-pointer*/
		arguments->Matrix[i] = allocateMemory(rows * sizeof(double*));

		for (j = 0; j < rows; j++)
		{
		    /*Speicherallokierung für die "rows" Zeilen*/
			arguments->Matrix[i][j] = arguments->M + (i * rows * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ********************************************************************************* */
/* initMatrices: Initialize chunk of matrix/matrices and some global variables       */
/* ********************************************************************************* */
static
void
initMatrices_parallel (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	uint64_t rows = arguments->rows;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;
	int offset = arguments->offset;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < rows; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0)
    {
        /*Nur der letzte und 0te Prozess müssen ihre erste (bzw letzte) Zeile initialisieren, aber alle Prozesse müssen den linken und rechten Rand initialisieren*/
        if(arguments->rang == 0)
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /*linker und rechter Rand*/
                    Matrix[g][i][0] = 1.0 - (h * i); //offset nicht nötig, da Rang = 0
                    Matrix[g][i][N] = h * i;
                }
                Matrix[g][0][N] = 0.0;

                for(i = 0; i <= N; i++)
                {
                    /*oberste Zeile*/
                    Matrix[g][0][i] = 1.0 - (h * i);
                }

            }
        }
        else if(arguments->rang == arguments->comm_size - 1)
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /*linker und rechter Rand*/
                    Matrix[g][i][0] = 1.0 - (h * (i + offset));
                    Matrix[g][i][N] = h * (i + offset);
                }

                for(i = 0; i <= N; i++)
                {
                    /*unterste Zeile*/
                    Matrix[g][rows - 1][i] = h * i;
                }
            }
        }
        else
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /*linker und rechter Rand*/
                    Matrix[g][i][0] = 1.0 - (h * (i + offset));
                    Matrix[g][i][N] = h * (i + offset);
                }
            }
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
    uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/* calculate_parallel: solves the equation in parallel                      */
/* ************************************************************************ */
static
void
calculate_parallel_jacobi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	int offset = arguments->offset;
	uint64_t rang = arguments->rang;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/*request-Variablen für das nichtblockierende Versenden von Nachrichten*/
    MPI_Request request_f; //first row
    MPI_Request request_l; //last row
    MPI_Status status_sendf;
    MPI_Status status_sendl;

    /*status-Variable zum blockierenden Empfangen von Nachrichten*/
    MPI_Status status_recv;

#ifdef HYBRID
    omp_set_dynamic(0);     //dynamische Teams ausschalten, damit runtime environment keinen Einfluss auf Anzahl Threads nimmt
    omp_set_num_threads(options->number); //setze Anzahl Threads auf vorgegebene Größe
#endif // HYBRID

    /* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;
#ifdef HYBRID
		/*
		Die durch openmp erstellten Threads teilen den i-Index untereinander auf, also die Zeilen der Matrix.
		Wir wÃ¤hlen den schedule-Modus guided, um eine bessere Lastverteilung zu erreichen. Die (hier: minimale) chunk-GrÃ¶ÃŸe von 1 ist normalerweise der default-Wert, wir setzen ihn aber der VollstÃ¤ndigkeit halber.
		Die Variablen i, j, star, residuum und maxResiduum werden von jedem Thread berechnet/verÃ¤ndert, weshalb wir sie auf private setzen.
		Die Ausnahme hierbei bildet maxResiduum. Dieses soll das Maximum aller Threads sein, weshalb wir nach Beendung des parallelen Blocks das Maximum Ã¼ber die reduction(max:)-Klausel berechnen.
		Die reduction-Klausel setzt auÃŸerdem die Variable maxResiduum implizit auf privat, weshalb wir es nicht in der private-Klausel haben.
		Alle anderen Variablen sind per default als shared gesetzt, weshalb wir die Klausel hier nicht explizit verwenden mÃ¼sen.
		*/
		#pragma omp parallel for schedule(guided, 1) private(i, j, star, residuum) reduction(max:maxResiduum)
#endif // HYBRID*
		/* over all rows */
		/* Der Index (rows - 1) entspricht der ersten Zeile des benachbarten Prozesses und soll nicht von diesem Thread bearbeitet werden.*/
		for (i = 1; (uint64_t)i < arguments->rows - 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + offset));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;

        if(rang == 0)
        {
            /*Prozess sendet seine letzte bearbeitete Zeile (nichtblockierend) an den nächsthöheren Prozess*/
            MPI_Isend((void*) &Matrix_Out[arguments->rows - 2][0], N + 1, MPI_DOUBLE, rang + 1, rang + 1, MPI_COMM_WORLD, &request_f);

            /*Prozess empfängt seine letzte Zeile (blockierend) vom nächsthöheren Prozess*/
            MPI_Recv((void*) &Matrix_Out[arguments->rows - 1][0], N + 1, MPI_DOUBLE, rang + 1, rang, MPI_COMM_WORLD, &status_recv);

            /*Prozess wartet auf das erfolgreiche Versenden seiner letzten bearbeiteten Zeile*/
            MPI_Wait(&request_f, &status_sendl);
        }
        else if(rang == arguments->comm_size - 1)
        {
            /*Prozess  sendet seine erste bearbeitete Zeile (nichtblockierend) an den vorherigen Prozess*/
            MPI_Isend((void*) &Matrix_Out[1][0], N + 1, MPI_DOUBLE, rang - 1, rang - 1, MPI_COMM_WORLD, &request_l);

            /*Prozess empfängt seine letzte Zeile (blockierend) vom vorherigen Prozess*/
            MPI_Recv((void*) &Matrix_Out[0][0], N + 1, MPI_DOUBLE, rang - 1, rang, MPI_COMM_WORLD, &status_recv);

            /*Prozess wartet auf das erfolgreiche Versenden seiner letzten bearbeiteten Zeile*/
            MPI_Wait(&request_l, &status_sendl);
        }
        else
        {
            /*Prozess sendet seine letzte bearbeitete Zeile (nichtblockierend) an den nächsthöheren Prozess*/
            MPI_Isend((void*) &Matrix_Out[arguments->rows - 2][0], N + 1, MPI_DOUBLE, rang + 1, rang + 1, MPI_COMM_WORLD, &request_l);
            /*Prozess  sendet seine erste bearbeitete Zeile (nichtblockierend) an den vorherigen Prozess*/
            MPI_Isend((void*) &Matrix_Out[1][0], N + 1, MPI_DOUBLE, rang - 1,  rang - 1, MPI_COMM_WORLD, &request_f);

            /*Prozess empfängt seine letzte Zeile (blockierend) vom vorherigen Prozess*/
            MPI_Recv((void*) &Matrix_Out[0][0], N + 1, MPI_DOUBLE, rang - 1, rang, MPI_COMM_WORLD, &status_recv);
            /*Prozess empfängt seine letzte Zeile (blockierend) vom nächsthöheren Prozess*/
            MPI_Recv((void*) &Matrix_Out[arguments->rows - 1][0], N + 1, MPI_DOUBLE, rang + 1, rang, MPI_COMM_WORLD, &status_recv);

            /*Prozess wartet auf das erfolgreiche Versenden*/
            MPI_Wait(&request_l, &status_sendl);
            MPI_Wait(&request_f, &status_sendf);
        }

        results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
		    /*Fehlerwert, welcher für die Abbruchbedingung genutzt wird.*/
            double lowest_precision = 100;

            /*Alle Prozesse senden synchron ihr maxResiduum and Prozess 0. Dieser speichert das größte residuum, also den größten Fehlerwert.*/
            MPI_Reduce((void*) &maxResiduum, &lowest_precision, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            /*Alle Prozesse empfangen synchron den größten Fehlerwert von Prozess 0.*/
            MPI_Bcast((void*) &lowest_precision, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (lowest_precision < options->term_precision)
			{
                term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/* calculate_parallel: solves the equation in parallel                      */
/* ************************************************************************ */
static
void
calculate_parallel_gauss (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;                   /* local variables for loops */
	int m1, m2;                 /* used as indices for old and new matrices */
	double star;                /* four times center value minus 4 neigh.b values */
	double residuum;            /* residuum of current iteration */
	double maxResiduum = 0;     /* maximum residuum value of a slave in iteration */
	double global_res = 0;      /* enthält das maximale Residuum der vorherigen Prozesse*/
	double highest_res = 0;     /* zum Versenden des maximalen Residuums*/
	int it_count = 0;           /* enthält die jetztige Iteration*/
	int flag_precision = 0;     /* Flag zur Überprüfung der Abbruchbedingung*/
	int last_it = 0;            /* enthält die Iteration in der abgebrochen wird*/
	int parity = 0;             /* für die Berechnung verschiedener Tags, um Iteration i und i+1 auseinander zu halten*/
	/*Da wir Isend verwenden um unsere erste berechnete Zeile zu versenden, kann es theoretisch dazu kommen dass 2 Nachrichten mit dem gleichen
	Tag versendet wurden, bevor eine ausgelesen wird. Um diese beiden Nachrichten auseinander halten zu können, wechseln wir mittels parity zwischen den
	Tags 12 und 13.*/

	int const N = arguments->N;
	double const h = arguments->h;

	int offset = arguments->offset;
	uint64_t rang = arguments->rang;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/*request-Variablen für das nichtblockierende Versenden/Empfangen von Nachrichten:
                        ..._f   first line
                        ..._l   last line
                        ..._r   residuum
                        ..._s   stop*/
    MPI_Request request_send_f;
    MPI_Request request_send_l;
    MPI_Request request_send_r;
    MPI_Request request_recv_l;
    MPI_Request request_recv_s;
    MPI_Status status_send_f;
    MPI_Status status_send_l;
    MPI_Status status_recv_s;

    /*status-Variable zum blockierenden Empfangen von Nachrichten*/
    MPI_Status status_recv;

	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;

		/*  EMPFANGEN DER ABBRUCHBEDINGUNG mit Irecv
        Falls Abbruch nach Iteration gewählt wurde, erwarten alle Prozesse den Wert der letzten zu berechnenden Iteration vom letzten Prozess.
        Dieser Wert wird jede Iteration geprüft, da wir den Abbruch mit Send verschicken (und ansonsten beim ersetzen von Send durch Ssend in einen
        Deadlock geraten).
        Die Nachricht kann erst empfangen werden, wenn der letzte Prozess angefangen hat zu rechnen (it_count + rang > arguments->comm_size - 1).
        Da die Abbruchbedingung nicht sofort benötigt wird, starten wir den Empfang mit Irecv.
		*/
		if((options->termination == TERM_PREC) && (rang != arguments->comm_size - 1) && (it_count + rang > arguments->comm_size - 1))
        {
            MPI_Irecv((void*) &last_it, 1, MPI_INT, arguments->comm_size - 1, 3, MPI_COMM_WORLD, &request_recv_s);
        }

		/*  EMPFANGEN LETZER ZEILE mit Irecv
		Alle Prozesse außer dem letzten empfangen die erste Zeile des nächsthöheren Prozesses und speichern sie als ihre letzte Zeile ab
		Da bei der 0. Iteration noch keine Kommunikation stattgefunden hat (it_count != 0), kann in diesem Fall auch nichts empfangen werden.
		Wir verwenden Irecv (nichtblockierend), da die letzte Zeile erst benötigt wird, wenn wir unsere vorletzte Zeile berechnen.
		Den Iterationscounter als Tag zu verwenden führt zu einem Fehler, da für den MPI-Tag nur 16 bits reserviert sind. Bei Iteration 131072 ist das Programm
		bei unseren Tests abgestürzt. Wir verwenden stattdessen parity, um Nachrichten aus verschiedenen Iterationen auseinander zu halten.
		Da wir die Nachricht aus der vorherigen Iteration erwarten, müssen wir für den Empfang den Wert von parity auf den Wert der vorherigen Iteration setzen,
		danach allerdings den Wert unserer jetztigen Iteration wiederherstellen.*/
		if((it_count != 0) && (rang != arguments->comm_size - 1))
        {
            parity = (parity + 1) % 2;
            MPI_Irecv((void*) &Matrix_In[arguments->rows - 1][0], N, MPI_DOUBLE, rang + 1, 12 + parity, MPI_COMM_WORLD, &request_recv_l);
            parity = (parity + 1) % 2;
        }

        /*  EMPFANGEN DES GLOBAL_RES    mit Recv
            EMPFANGEN DER 0. ZEILE      mit Recv
        Alle Prozesse außer dem ersten empfangen die letzte Zeile des vorherigen Prozesses und speichern sie als ihre neue erste Zeile ab.
        Vor dem Empfang kann die Matrix nicht korrekt berechnet werden, weshalb wir hier Recv (blockierend) verwenden.
        Des weiteren empfangen wir das größte Residuum aller vorherigen Prozesse.*/
		if(rang != 0)
        {
            if(options->termination == TERM_PREC)
            {
                MPI_Recv((void*) &global_res, 1, MPI_DOUBLE, rang - 1, 11, MPI_COMM_WORLD, &status_recv);
            }
            MPI_Recv((void*) &Matrix_In[0][0], N, MPI_DOUBLE, rang - 1, 10, MPI_COMM_WORLD, &status_recv);
        }

        for (i = 1; (uint64_t)i < arguments->rows - 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + offset));
			}

			/*  SENDEN DER ZWEITEN ZEILE MIT ROWS >= 4 mit Isend
			Alle Prozesse außer dem ersten senden ihre Zeile 1 an den vorherigen Prozess.
			Wir verwenden Irecv (nichtblockierend), da es für diesen Prozess egal ist wann genau das Senden stattfindet.
			Wir müssen allerdings beachten, dass wir die Zeile fertig berechnet haben. Wenn (i == 2) gilt, ist die Berechnung der ersten Zeile abgeschlossen.
			Wir verwenden Isend, da es für diesen Prozess (fast) egal ist, wann genau das Versenden stattfindet.
			Wir verwenden parity, um beim Empfang verschiedene Iterationen auseinander halten zu können.*/
			if((rang != 0) && (i == 2))
            {
                MPI_Isend((void*) &Matrix_Out[1][0], N, MPI_DOUBLE, rang - 1, 12 + parity, MPI_COMM_WORLD, &request_send_f);
            }

            /*  ABSCHLUSS SENDEN DER VORLETZTEN ZEILE mit Isend
                ABSCHLUSS EMPFANGEN DER LETZTEN ZEILE mit Irecv
            Alle Prozesse außer dem letzten schließen das Senden ihrer vorletzten Zeile und das Empfangen ihrer letzten Zeile ab.
            In Iteration 0 wurde noch nicht versendet, weshalb in diesem Fall nicht gewartet werden muss.
            Falls (i == arguments->rows - 2) gilt verändern wir als nächstes die vorletzte Zeile, weshalb sowohl das Versenden als auch der Empfang
            vorher abgeschlossen werden müssen.*/
            if((it_count != 0) && (rang != arguments->comm_size - 1) && ((uint64_t)i == arguments->rows - 2))
            {
                MPI_Wait((void*) &request_send_l, &status_send_l);
                MPI_Wait((void*) &request_recv_l, &status_send_l);
            }

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		/*  SENDEN DER ZWEITEN ZEILE MIT ROWS < 4 mit Isend
        Alle Prozesse außer dem ersten senden ihre Zeile 1 an den vorherigen Prozess.
		Wir verwenden Irecv (nichtblockierend), da es für diesen Prozess egal ist wann genau das Senden stattfindet.
		Diese if-Abfrage ist für das Abfangen von Randfällen gedacht. Bei ausreichend großen Matrizen wird die Zeile 1
		bereits in der for-Schleife versendet. Falls die Matrix nicht groß genug ist (rows < 4) soll dies hier nachgeholt
		werden. Wir können also die selben Parameter (request_send_f, tag) verwenden, ohne in einen Konflikt zu geraten.*/
        if((rang != 0) && (arguments->rows < 4))
        {
            MPI_Isend((void*) &Matrix_Out[1][0], N, MPI_DOUBLE, rang - 1, 12 + parity, MPI_COMM_WORLD, &request_send_f);
        }

		/*Alle Prozesse außer dem letzten schließen das Versenden ihres highest_res ab, bevor dieses verglichen wird.
		In Iteration 0 hat noch keine Kommunikation stattgefunden, weshalb hier nicht gewartet werden muss (it_count !=0 ).*/
        if((options->termination == TERM_PREC) && (it_count != 0) && (rang != arguments->comm_size - 1))
        {
            MPI_Wait(&request_send_r, &status_send_l);
        }

		/* Alle Prozesse gleichen ihr lokales maxResiduum mit dem des vorherigen Prozesses ab, und speichern das größere.*/
		if(options->termination == TERM_PREC)
        {
            highest_res = (maxResiduum > global_res) ? maxResiduum : global_res;
        }

        /*  SENDEN DER VORLETZTEN ZEILE mit Isend
            SENDEN DES MAXIMALEN RESIDUUMS mit Isend
        Alle Prozesse außer dem letzten versenden ihre vorletzte Zeile an den nächsten Prozess.
        Alle Prozesse außer dem letzten versenden das maximale Residuum an den nächsten Prozess.
        (entweder ihr eigenes, oder das Residuum das ihnen zugesendet wurde)
        Wann genau versendet wird ist diesem Prozess egal, weshalb wir Isend verwenden.*/
        if(rang != arguments->comm_size - 1)
        {
            if(options->termination == TERM_PREC)
            {
                MPI_Isend((void*) &highest_res, 1, MPI_DOUBLE, rang + 1, 11, MPI_COMM_WORLD, &request_send_r);
            }
            MPI_Isend((void*) &Matrix_Out[arguments->rows - 2][0], N, MPI_DOUBLE, rang + 1, 10, MPI_COMM_WORLD, &request_send_l);
        }

        /*  ABSCHLUSS: SENDEN DER ZWEITEN ZEILE mit Isend
		Alle Prozesse außer dem ersten senden ihre erste berechnete Zeile an den vorherigen Prozess.
		Das Senden muss hier abgeschlossen werden, da wir diese Zeile beim Start der nächsten Iteration verändern.
		Wir hätten das Senden auch am Anfang der while-Schleife abschließen können, allerdings hätte dies nicht zur
		besseren Verständlichkeit des Codes beigetragen, da wir für den Abschluss in der letzten Iteration noch ein
		weiteres Wait hätten benutzen müssen.*/
        if(rang != 0)
        {
            MPI_Wait(&request_send_f, &status_send_f);
        }

        /*  ABSCHLUSS: EMPFANGEN DER ABBRUCHBEDINGUNG mit Irecv
        Die Abbruchbedingung wir als nächstes geprüft, also muss die Nachricht bis hier hin empfangen worden sein.*/
        if((options->termination == TERM_PREC) && (rang != arguments->comm_size - 1) && (it_count + rang > arguments->comm_size - 1))
        {
            MPI_Wait(&request_recv_s, &status_recv_s);
        }

		/*Iterationscounter pro Iteration um einen inkrementieren*/
		it_count++;

		parity = (parity + 1) % 2;

		results->stat_precision = maxResiduum;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
		    /*Der letzte Prozess überprüft ob die Abbruchbedingung erreicht wurde.
		    Da alle vorherigen Prozesse das jeweils größte maxResiduum weiterreichen, erhält der letzte
		    Prozess immer das größte maxResiduum der ganzen Matrix.
		    Falls Abbruchbedingung erreicht wurde, sendet er eine Nachricht an Prozess 0. Deren Inhalt ist dabei egal.*/
		    if(rang == arguments->comm_size - 1)
            {
                if((!flag_precision) && (global_res < options->term_precision))
                {
                    last_it = it_count + arguments->comm_size;
                    flag_precision = 1;

                }

                for(i = 0; i < (int)arguments->comm_size - 1; i++)
                {
                    if((last_it == 0) || (it_count + (int)(arguments->comm_size - 1) - i < last_it))
                    {
                            MPI_Send((void*) &last_it, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
                    }
                }
            }

            /*Falls wahr wird die while-Schleife verlassen.*/
            if(it_count == last_it)
            {
                term_iteration = 0;
            }
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	/*  ABSCHLUSS: SENDEN DER VORLETZTEN ZEILE mit Isend
	    ABSCHLUSS: SENDEN DES MAXIMALEN RESIDUUMS mit Isend
	Für die letzte Iteration muss das letzte Versenden der vorletzen Zeile und des Residuums abgeschlossen werden. Dieses Versenden wird normalerweise
	innerhalb der äußeren for-Schleife abgeschlossen, welche wir nach der letzten Iteration allerdings nicht mehr betreten.*/
	if(rang != arguments->comm_size - 1)
    {
        MPI_Wait(&request_send_l, &status_send_l);

        if(options->termination == TERM_PREC)
        {
            MPI_Wait(&request_send_r, &status_send_l);
        }
    }

    /*Debugging und testen*/
	//printf("rang: %" PRIu64 " beendet\niteration:       %d\nlokales residuum:     %f\nglobales residuum:    %f\nterm_precision:       %f\n\n", rang, it_count, maxResiduum, global_res, options->term_precision);
    results->stat_iteration = it_count;

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/* ********************************************************************************* */
/*  displayStatistics_parallel: displays some statistics about the calculation       */
/* ********************************************************************************* */
static
void
displayStatistics_parallel (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	double speicher = (N + 1) * arguments->rows * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0;
	double speicherGesamt;

	/*Alle Prozesse senden synchron ihren Speicherbedarf an Prozess 0. Dieser speichert die Summe der Werte.*/
    MPI_Reduce((void*) &speicher, &speicherGesamt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double fehler = 0.0;

    /*Alle Prozesse senden synchron ihren Fehler an Prozess 0. Dieser speichert das Maximum.*/
    MPI_Reduce((void*) &results->stat_precision, &fehler, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if(arguments->rang == 0)
    {
        printf("Berechnungszeit:    %f s \n", time);
        printf("Speicherbedarf:     %f MiB\n", speicherGesamt);
        printf("Berechnungsmethode: ");

        if (options->method == METH_GAUSS_SEIDEL)
        {
            printf("Gauß-Seidel");
        }
        else if (options->method == METH_JACOBI)
        {
            printf("Jacobi");
        }

        printf("\n");
        printf("Interlines:         %" PRIu64 "\n",options->interlines);
        printf("Stoerfunktion:      ");

        if (options->inf_func == FUNC_F0)
        {
            printf("f(x,y) = 0");
        }
        else if (options->inf_func == FUNC_FPISIN)
        {
            printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
        }

        printf("\n");
        printf("Terminierung:       ");

        if (options->termination == TERM_PREC)
        {
            printf("Hinreichende Genaugkeit");
        }
        else if (options->termination == TERM_ITER)
        {
            printf("Anzahl der Iterationen");
        }

        printf("\n");
        printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
        printf("Norm des Fehlers:   %e\n", fehler);
        printf("\n");
    }
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix_parallel                       **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix_parallel (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
    int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	MPI_Status status;

	/*Es sollen 9 Zeilen ausgegeben werden. Hierfür startet Prozess 0 mit y = 0 und verlässt die for-Schleife, sobald
	(y * (interlines + 1)) außerhalb seines Bereiches liegt. Dann wird der aktuelle y-Wert an den nächsten Prozess gesendet.
	Jeder Prozess überprüft, ob (y * (interlines + 1)) in dem Bereich der von ihm berechneten Zeilen liegt, gibt ggf die betreffende
	Zeile aus und sendet den aktuellen y-Wert an den nächsthöheren Prozess.*/
	if(arguments->rang == 0)
    {
        for (y = 0; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->rows - 1; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1))][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);

        MPI_Send((void*) &y, 1, MPI_INT, arguments->rang + 1, 6, MPI_COMM_WORLD);
    }
    else if(arguments->rang == arguments->comm_size - 1)
    {
        MPI_Recv((void*) &y, 1, MPI_INT, arguments->rang - 1, 6, MPI_COMM_WORLD, &status);

        for (; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->offset + arguments->rows && (y * (interlines + 1)) > arguments->offset; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) - arguments->offset][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);
    }
    else
    {
        MPI_Recv((void*) &y, 1, MPI_INT, arguments->rang - 1, 6, MPI_COMM_WORLD, &status);

        for (; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->offset + arguments->rows - 1 && (y * (interlines + 1)) > arguments->offset; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) - arguments->offset][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);

        MPI_Send((void*) &y, 1, MPI_INT, arguments->rang + 1, 6, MPI_COMM_WORLD);
    }
}

static
void
sequential_block(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int argc, char** argv)
{
    askParams(options, argc, argv);
    initVariables(arguments, results, options);

    allocateMatrices(arguments);
	initMatrices(arguments, options);

	gettimeofday(&start_time, NULL);
	calculate(arguments, results, options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(arguments, results, options);
	displayMatrix(arguments, results, options);

	freeMatrices(arguments);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
#ifdef PARALLEL
	MPI_Init(NULL, NULL);

    /*Abfragen und Abspeicherung des Prozessrangs*/
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    arguments.rang = world_rank;

    /*Abfrage und Speicherung der Anzahl der Prozesse insgesamt*/
	int world_size = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	arguments.comm_size = world_size;

	if(world_size > 1)
    {
        /*Prozess 0 nimmt die Parameter entgegen*/
        if(world_rank == 0)
        {
            askParams(&options, argc, argv);
        }

        /*Prozess 0 übermittlelt die Parameter an alle anderen Prozesse*/
        MPI_Bcast((void*) &options.number, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.method, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.interlines, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.inf_func, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.termination, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.term_iteration, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*) &options.term_precision, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        initVariables(&arguments, &results, &options);

        /*Berechnen der Anzahl der Zeilen, die je ein Prozess bearbeitet.
        Der Rückgabewert ist 0 falls die Aufteilung erfolgreich war, andernfalls wird das Programm beendet.*/
        if(calc_chunk(&arguments) != 0)
        {
            MPI_Finalize();

            return 0;
        }

        allocateMatrices_parallel(&arguments);

        initMatrices_parallel(&arguments, &options);

        /*Synchronisation der Prozesse für einheitliche Messungen.*/
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&start_time, NULL);

        if (options.method != METH_JACOBI)
        {
            calculate_parallel_gauss(&arguments, &results, &options);
        }
        else
        {
            calculate_parallel_jacobi(&arguments, &results, &options);
        }

        /*Synchronisation der Prozesse für einheitliche Messungen.*/
        MPI_Barrier(MPI_COMM_WORLD);
        gettimeofday(&comp_time, NULL);

        displayStatistics_parallel(&arguments, &results, &options);
        displayMatrix_parallel(&arguments, &results, &options);

        freeMatrices(&arguments);
    }
    else
    {
        sequential_block(&arguments, &results, &options, argc, argv);
    }

	MPI_Finalize();
#else
    sequential_block(&arguments, &results, &options, argc, argv);
#endif
	return 0;
}
