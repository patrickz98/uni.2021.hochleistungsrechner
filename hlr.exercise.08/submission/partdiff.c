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
	int offset;
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
void
calc_chunk (struct calculation_arguments* arguments)
{
    uint64_t const N = arguments->N;

    int world_size = arguments->comm_size;

	/*Berechnung des chunks den je ein Prozess übernimmt. Entspricht der Anzahl an Zeilen der Matrix die vom jeweiligen Prozess gehalten werden, NICHT der Anzahl der zu berechnenden Zeilen.*/
	uint64_t rows = (N + 1) / world_size;

	/*Da sich die Matrix nur im Ausnahmefall glatt auf alle Prozesse aufteilen lässt, müssen einige Prozesse mehr Zeilen
	übernehmen als andere. In diesem Fall übernehmen die ersten x Prozesse eine Zeile mehr.*/
	if((N + 1) % world_size != 0)
    {
        if(arguments->rang < ((N + 1) % world_size))
        {
            rows++;
        }
    }

    /*Da jeder Prozess die unterste Zeile des vorherigen Prozesses und die erste Zeile des nächsten Prozesses für die Berechnung braucht,
	werden zusätlich zum chunk 2 weitere Zeilen allokiert. Der erste und letzte Prozess benötigen nur eine zusätzliche Zeile.*/
	if(arguments->rang == 0 || (arguments->rang + 1) == arguments->comm_size)
    {
        rows++;
    }
    else
    {
        rows += 2;
    }

    /*offset Beschreibt die Position der jeweiligen Zeile in der Gesamtmatrix. Für die sin-Berechnung bei clalculate_parallel benötigt,
    sowie bei init_matrices_parallel und display_matrices_parallel.
    Hierfür müssen wir mit 1 subtrahieren, da die letzte Zeile des n-ten Prozesses der ersten Zeile des (n-1)-ten Prozesses entspricht.
    Da sich die Zeilen meistens nicht glatt aufteilen lassen, müssen wir die übrig gebliebenen Zeilen dazuaddieren.*/
	if(arguments->rang == 0)
    {
        arguments->offset = 0;
    }
    else
    {
        arguments->offset = rows * arguments->rang - 1 + (arguments->rang <= (N + 1) % world_size ? arguments->rang : (N + 1) % world_size);
    }

    arguments->rows = rows;
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
        else if(arguments->rang == arguments->comm_size)
        {
            for(g = 0; g < arguments->num_matrices; g++)
            {
                for(i = 0; i < rows; i++)
                {
                    /*linker und rechter Rand*/
                    Matrix[g][i][0] = 1.0 - (h * (i + offset));
                    Matrix[g][i][N] = h * (i + offset);
                }
                Matrix[g][rows][0] = 0.0;

                for(i = 0; i <= N; i++)
                {
                    /*unterste Zeile*/
                    Matrix[g][rows][i] = h * i;
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
calculate_parallel (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
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

        /*Fehlerwert, welcher für die Abbruchbedingung genutzt wird.*/
        double lowest_precision = 100;

        /*Alle Prozesse senden synchron ihr maxResiduum and Prozess 0. Dieser speichert das größte residuum, also den größten Fehlerwert.*/
        MPI_Reduce((void*) &maxResiduum, &lowest_precision, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        /*Alle Prozesse empfangen synchron den größten Fehlerwert von Prozess 0.*/
        MPI_Bcast((void*) &lowest_precision, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        results->stat_precision = lowest_precision;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
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

	/*Alle Prozesse senden synchron ihren Speicherbedarf and Prozess 0. Dieser speichert die Summe der Werte.*/
    MPI_Reduce((void*) &speicher, &speicherGesamt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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
        printf("Norm des Fehlers:   %e\n", results->stat_precision);
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

	if(arguments->rang == 0)
    {
        for (y = 0; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->rows; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) % arguments->rows][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);

        MPI_Ssend((void*) &y, 1, MPI_INT, arguments->rang + 1, 3, MPI_COMM_WORLD);
    }
    else if(arguments->rang == arguments->comm_size - 1)
    {
        MPI_Recv((void*) &y, 1, MPI_INT, arguments->rang - 1, 3, MPI_COMM_WORLD, &status);

        for (; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->offset + arguments->rows; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) % arguments->rows][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);
    }
    else
    {
        MPI_Recv((void*) &y, 1, MPI_INT, arguments->rang - 1, 3, MPI_COMM_WORLD, &status);

        for (; y < 9 && (uint64_t)(y * (interlines + 1)) < arguments->offset + arguments->rows; y++)
        {
            for (x = 0; x < 9; x++)
            {
                printf ("%7.4f", Matrix[(y * (interlines + 1)) % arguments->rows][x * (interlines + 1)]);
            }

            printf ("\n");
        }

        fflush (stdout);

        MPI_Ssend((void*) &y, 1, MPI_INT, arguments->rang + 1, 3, MPI_COMM_WORLD);
    }
}

static
void
sequential_block(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
    allocateMatrices(arguments);
	initMatrices(arguments, options);

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

        /*Falls Gauss-Seidel-Verfahren gewählt wurde, führt Prozess 0 dieses aus. Die restlichen Prozesse warten bis gemeinsam beendet wird.*/
        if (options.method != METH_JACOBI)
        {
            if(world_rank == 0)
            {
                sequential_block(&arguments, &results, &options);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Finalize();

            return 0;
        }

        initVariables(&arguments, &results, &options);

        /*Berechnen der Anzahl der Zeilen, die je ein Prozess bearbeitet*/
        calc_chunk(&arguments);

        allocateMatrices_parallel(&arguments);

        initMatrices_parallel(&arguments, &options);

        gettimeofday(&start_time, NULL);

        calculate_parallel(&arguments, &results, &options);

        /*Barriere um Messung einheitlich zu halten*/
        MPI_Barrier(MPI_COMM_WORLD);

        gettimeofday(&comp_time, NULL);

        displayStatistics_parallel(&arguments, &results, &options);
        displayMatrix_parallel(&arguments, &results, &options);

        freeMatrices(&arguments);
    }
    else
    {
        askParams(&options, argc, argv);
        initVariables(&arguments, &results, &options);
        sequential_block(&arguments, &results, &options);
    }

	MPI_Finalize();
#else
    askParams(&options, argc, argv);
    initVariables(&arguments, &results, &options);
    sequential_block(&arguments, &results, &options);
#endif
	return 0;
}
