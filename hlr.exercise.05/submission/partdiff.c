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
#include <pthread.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

typedef enum {false, true} bool;

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

// Mutex für parallele Verarbeitung
pthread_mutex_t mutex_maxResiduum;

// maxResiduum ist jetzt global, damit alle Threads darauf zugreifen können
double maxResiduum;

/*struct um Threads Variablen zu übergeben*/
struct thread_variables
{
    int id;
    int start_index, end_index; //Indizes für äußere loop
    int N; //Anzahl Elemente einer Zeile der Matrix, für innere loop
    double** Matrix_In;
    double** Matrix_Out;
    bool res; //um if-Abfragen zu beschleunigen
    bool func; //um if-Abfragen zu beschleunigen
    double pih;
	double fpisin;
};

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

static
void*
for_loops(void* t_var_array)
{
    struct thread_variables* vars = (struct thread_variables*) t_var_array;

    double residuum = 0.0;
    double mResiduum = 0.0;
    double star;

    /* over all rows */
    for (int i = vars->start_index; i < vars->end_index; i++)
	{
		// fprintf(stderr, "######### [%d] i = %d\n", vars->id, i);

		double fpisin_i = 0.0;

		if (vars->func)
		{
			fpisin_i = vars->fpisin * sin(vars->pih * (double)i);
		}

		/* over all columns */
		for (int j = 1; j < vars->N; j++)
		{
			star = 0.25 * (vars->Matrix_In[i-1][j]
				+ vars->Matrix_In[i][j-1]
				+ vars->Matrix_In[i][j+1]
				+ vars->Matrix_In[i+1][j]);

			if (vars->func)
			{
				star += fpisin_i * sin(vars->pih * (double)j);
			}

			if (vars->res)
			{
				residuum = vars->Matrix_In[i][j] - star;
				residuum = (residuum < 0) ? -residuum : residuum;
                mResiduum = (residuum < mResiduum) ? mResiduum : residuum;
			}

			vars->Matrix_Out[i][j] = star;
		}
	}

	if (vars->res)
    {
        // Versuche lock auf maxResiduum zu erhalten
        pthread_mutex_lock(&mutex_maxResiduum);

        // Verändere maxResiduum
        maxResiduum = (mResiduum < maxResiduum) ? maxResiduum : mResiduum;

        // Gebe lock auf maxResiduum frei
        pthread_mutex_unlock(&mutex_maxResiduum);
    }

	pthread_exit((void*) t_var_array);
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_parallel (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	// used as indices for old and new matrices
	int m1, m2;

	// Für die Ausgabe des Fehlercodes, sollte einer auftreten. Für debugging
	int rc;

	// Array mit allen Threads
    pthread_t my_threads[options->number];

	// Wir kreieren die Threads explizit als joinable,
	// um Fehler aufgrund Umgebungsvariablen zu vermeiden
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    /*Wir initialisieren die Mutex-Variable*/
    pthread_mutex_init(&mutex_maxResiduum, NULL);

	int const N = arguments->N;
	double const h = arguments->h;

	//Wir speichern das Ergebnis als boolean, um in den verschachtelten for-loops nicht immer die gleiche Abfrage starten zu müssen
	bool func = options->inf_func == FUNC_FPISIN;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	int thread_limit;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
		thread_limit = options->number;
	}
	else
	{
		m1 = 0;
		m2 = 0;
		thread_limit = 1;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	// Für die Berechnung der zu bearbeitenden Indizes der jeweiligen Threads
	int chunk = N / thread_limit;

	// For-loop counter variable
	int inx = 0;

	// This thread is needed to store thread variables.
    struct thread_variables t_var_array[thread_limit];

	while (term_iteration > 0)
	{
		/*if (term_iteration % 100 == 0)
		{
			fprintf(stderr, "######### iteration: %d\n", term_iteration);
		}*/

		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0.0;

		// Wir speichern das Ergebnis als boolean,
		// um in den verschachtelten for-loops nicht
		// immer die gleiche Abfrage starten zu müssen
		bool res = false;
		if (options->termination == TERM_PREC || term_iteration == 1)
        {
            res = true;
        }

		// fprintf(stderr, "######### N: %d\n", N);
		// fprintf(stderr, "######### chunk: %d\n", chunk);
		// fprintf(stderr, "######### threads: %d\n", options->number);

		// Create Threads
		for (inx = 0; inx < thread_limit; inx++)
        {
            //Der letzte Thread muss etwas mehr arbeiten, falls sich die Matrix nicht perfekt aufteilen lässt
			unsigned int end_index = (inx == (thread_limit - 1))
				? (unsigned int) N
				: (unsigned int) ((inx + 1) * chunk) + 1;

			// Initialisierung der Variablen
			struct thread_variables vars = {
				.N           = N,
				.start_index = (inx * chunk) + 1,
				.end_index   = end_index,
				.Matrix_In   = Matrix_In,
				.Matrix_Out  = Matrix_Out,
				.res         = res,
				.func        = func,
				.fpisin      = fpisin,
				.pih         = pih,
				.id          = inx,
			};

			t_var_array[inx] = vars;

			// fprintf(stderr, "######### [%d] start_index: %4d\n", inx, vars.start_index);
			// fprintf(stderr, "######### [%d] end_index:   %4d\n", inx, vars.end_index);

			rc = pthread_create(&my_threads[inx], &attr, for_loops, (void*) &t_var_array[inx]);
            if (rc)
            {
                fprintf(stderr, "ERROR: [%d] return code from pthread_create() is %d\n", inx, rc);
                exit(-1);
            }
        }

		// exit(0);

        // Wir führen die Threads zusammen, um sicher zu sein dass alle Threads fertig sind.
		// Andernfalls kann es zu Datenabhängigkeiten kommen
        for (inx = 0; inx < thread_limit; inx++)
        {
            rc = pthread_join(my_threads[inx], (void*) NULL);
            if (rc)
            {
                fprintf(stderr, "ERROR: [%d] return code from pthread_join() is %d\n", inx, rc);
                exit(-1);
            }
        }

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		int inx = m1;
		m1 = m2;
		m2 = inx;

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

	/*Wir geben das Attribut und die Mutex wieder frei*/
	pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&mutex_maxResiduum);
	results->m = m2;
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

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
#if defined(POSIX)
	calculate_parallel(&arguments, &results, &options);
#else
	calculate(&arguments, &results, &options);
#endif
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
