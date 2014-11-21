/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
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

#include "partdiff-posix.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
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
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


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
		printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", size);
		/* exit program */
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
	uint64_t g, i, j;                                /*  local variables for loops   */

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

//Struktur, um der Funktion do_calc Parameter zu übergeben
struct pthread_calc_args {
    //Der Rückgabewert der Funktion do_calc
    double returnValue;
    
    //Variablen für den Thread-Bereich
    int thread_number;
    int start_i;
    int end_i;
    
    //Parameter der Funktion calculate, welche zwecks der Verwendung in do_calc
    //ausgelagert wurden
    double** Matrix_Out;
    double** Matrix_In;
    int term_iteration;
    int inf_func;
    int termination;
    double pih;
    double fpisin;
    
    double h;
    int N;
};

/**
 * Funktion für pThread, die die Berechnung eines Datenblocks übernimmt.
 * @param t_args Erwartet struct pthread_calc_args* mit den notwendigen Daten
 * @return NULL
 */
void *do_calc(void *t_args)
{
    
        struct pthread_calc_args *args = (struct pthread_calc_args*) t_args;
	double local_maxresiduum = 0.0;
	int i, j;
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
        
        //printf("Hello from thread %i (data %i to %i)\n", args->thread_number, args->start_i, args->end_i);
    
        /* over all rows */
        //for (i = 1; i < args->N; i++)
        for(i = args->start_i; i <= args->end_i; i++)
        {
                double fpisin_i = 0.0;

                if (args->inf_func == FUNC_FPISIN)
                {
                        fpisin_i = args->fpisin * sin(args->pih * (double)i);
                }

                /* over all columns */
                for (j = 1; j < args->N; j++)
                {
                        star = 0.25 * (args->Matrix_In[i-1][j] + args->Matrix_In[i][j-1] + args->Matrix_In[i][j+1] + args->Matrix_In[i+1][j]);

                        if (args->inf_func == FUNC_FPISIN)
                        {
                                star += fpisin_i * sin(args->pih * (double)j);
                        }

                        if (args->termination == TERM_PREC || args->term_iteration == 1)
                        {
                                residuum = args->Matrix_In[i][j] - star;
                                residuum = (residuum < 0) ? -residuum : residuum;
                                local_maxresiduum = (residuum < local_maxresiduum) ? local_maxresiduum : residuum;
                        }

                        args->Matrix_Out[i][j] = star;
                }
        }
    
        //printf("Finished from thread %i (data %i to %i) (%ld)\n", args->thread_number, args->start_i, args->end_i, local_maxresiduum);
        
        //Rückgabewert speichern
        args->returnValue = local_maxresiduum;
    
        return NULL;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
	int i;//, j; /* local variables for loops */
        int m1, m2; /* used as indices for old and new matrices */
        //double star; /* four times center value minus 4 neigh.b values */
        //double residuum; /* residuum of current iteration */
        double maxresiduum; /* maximum residuum value of a slave in iteration */
        //int const N = arguments->N;
        double const h = arguments->h;
        double pih = 0.0;
        double fpisin = 0.0;
        int term_iteration = options->term_iteration;
        
	int rc;

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

        //////////////////////////////////////////////7
        /// Pthread
        // Jeweils Arrays, als Speichertypen
        
        //Struktur, um die einzelnen Thread Argumente zu speichern.
        struct pthread_calc_args thread_args[options->number]; 
        //Array um die pthread_t Daten zu speichern
	pthread_t threads[options->number];
        //Die 'Breite' eines Thread-Datenblocks bestimmen
        double data_width = (double)(arguments->N - 1) / options->number;
        
        
	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];
                
                maxresiduum = 0.0;

                //printf("ITER: %4i START\n", term_iteration);
                
                //Jeden Thread initalisieren und starten
                for(i=0; i<options->number; i++) {
                    //Die Threadargumente initalisieren
                    thread_args[i].thread_number = i;
        
                    thread_args[i].start_i = data_width * i + 1; //Erster Index des Threadsdatenblocks
                    thread_args[i].end_i = data_width * (i + 1); //Letzter Index des Threadsdatenblocks
                    
                    thread_args[i].Matrix_In = Matrix_In;
                    thread_args[i].Matrix_Out = Matrix_Out;
                    thread_args[i].term_iteration = term_iteration;
                    thread_args[i].inf_func = options->inf_func;
                    thread_args[i].termination = options->termination;
                    thread_args[i].pih = pih;
                    thread_args[i].fpisin = fpisin;
                    
                    thread_args[i].h = arguments->h;
                    thread_args[i].N = arguments->N;
                    
                    //Den Thread starten
                    rc = pthread_create(&threads[i], NULL, do_calc, &thread_args[i]);
                    if (rc)
                    {
                        printf("Error: Thread %i could not be created\n", i);
                    }
                }
                
                //Auf die einzelnen Threads warten
                for(i=0; i<options->number; i++) {
                    pthread_join(threads[i], NULL);
                    //Den 'Rückgabewert' maxresiduum auswerten und 'zusammenfassen'
                    maxresiduum = (maxresiduum < thread_args[i].returnValue) ? thread_args[i].returnValue : maxresiduum;
                }
                //printf("ITER: %4i END\n", term_iteration);

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
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
		printf("Gauss-Seidel");
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
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
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

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */

	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);                   /*  start timer         */
	calculate(&arguments, &results, &options);                                      /*  solve the equation  */
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);                                       /*  free memory     */

	return 0;
}
