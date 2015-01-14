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
#include <mpi.h>

#include "partdiff-par.h"
/*
INFO:

Isend and Irecv aren't used so far, because for easier debugging.
We tried first to use Isend and Irecv, but the code wasn't working that time.

All "INFO..."-macros are used for debugging only.
If DEBUG is not defined, these macros are doing nothing.

Notice ourself:
Try to improve the code and performance by using those.
*/

/*
TODOS:
 - check the termination after precision. Its getting close, but the parallel program is somehow getting to the result with ~10% less iterations.
*/

#define UNUSED(x) (void)x;

//#define DEBUG
#ifdef DEBUG
    #define INFO_STR(rank, size, x) printf("rank %2d of %2d (%3d): %s\n", rank+1, size, __LINE__, x)
    #define INFO_DOUBLE(rank, size, number, x, y) printf("rank %2d of %2d (%3d) (%2d|%2d): %f\n", rank+1, size, __LINE__, x, y, number)
	#define INFO_DOUBLESTR(rank, size, str, x) printf("rank %2d of %2d (%3d): " str "\n", rank+1, size, __LINE__, x)
    #include <string.h>
    #define INFO_DATA(rank, size, data, datasize) \
            do { \
                char* str = (char*)malloc(( (datasize+1) * 8 + 1) * sizeof(char)); \
                for(int _i_=0; _i_<datasize; _i_++) { \
                    sprintf(str + (8*_i_), ", %1.4f", data[_i_]); \
                } \
                INFO_STR(rank, size, str); \
    			free(str); \
            } while(0)
#else
    #define INFO_STR(rank, size, x) UNUSED(x)
    #define INFO_DOUBLE(rank, size, number, x, y) UNUSED(x)
	#define INFO_DOUBLESTR(rank, size, str, x) UNUSED(x)
    #define INFO_DATA(rank, size, data, datasize) UNUSED(data)
#endif

#define MPI_TAG_UP   333
#define MPI_TAG_DOWN 666
#define MPI_TAG_ABORT 777
#define MPI_TAG_PRECISION 999

struct calculation_arguments {
    uint64_t N;             /* number of spaces between lines (lines=N+1) in total */
    uint64_t num_matrices;  /* number of matrices                                  */
    double h;               /* length of a space between two lines                 */
    double ***Matrix;       /* index matrix used for addressing M                  */
    double *M;              /* two matrices with real values                       */

    uint64_t N_rank;        /* like N, but for this SPECIFIC Process N_rank <= rank */
    int row_start;          /* the first row/line for this rank                     */
    int row_end;            /* the last row/line for this rank                      */
};

struct calculation_results {
    uint64_t m;
    uint64_t stat_iteration;/* number of current iteration                 */
    double stat_precision;  /* actual precision of all slaves in iteration */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;  /* time when program started                      */
struct timeval comp_time;   /* time when calculation completed                */



/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables(struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options) {
    arguments->N = (options->interlines * 8) + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = 1.0 / arguments->N;

    results->m = 0;
    results->stat_iteration = 0;
    results->stat_precision = 0;

    //used for the distribution of the remaining lines to some ranks, in case N + 1 % size != 0
    int rest = (arguments->N + 1) % options->size;

    //calculates the number of lines for the specific rank
    arguments->N_rank = (arguments->N + 1) / options->size;

    //for the correction of the row_start/row_end variable
    //taking into account the rest_offset(s) of the previous ranks
    int rest_offset = 0;
    for (int i = 0; i < options->rank; ++i) {
        if (i < rest) {
            rest_offset++;
        }
    }
    // compute start rank + all offsets off the smaller ranks
    arguments->row_start = options->rank * arguments->N_rank + rest_offset;
    //offset of my rank
    if (options->rank < rest) {
        rest_offset++;
    }
    arguments->row_end = (options->rank + 1) * arguments->N_rank + rest_offset - 1;
    //give every process his first and his last line, that he doesnt compute
    //but needs from the other processes to compute his lines
    
    //for the first rank
    if (options->rank > 0) {
        arguments->row_start--;
        arguments->N_rank++;
    }
    //for the last rank
    if (options->rank < (options->size - 1)) {
        arguments->row_end++;
        arguments->N_rank++;
    }
    //distributes the remaining lines to the first ranks
    if (options->rank < rest) {
        arguments->N_rank++;
    }
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices(struct calculation_arguments* arguments) {
    uint64_t i;

    for (i = 0; i < arguments->num_matrices; i++) {
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
allocateMemory(size_t size) {
    void *p;

    if ((p = malloc(size)) == NULL) {
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
allocateMatrices(struct calculation_arguments* arguments) {
    /* M can be described as a three dimensional cube (Array)

           +----------+ 
          1.         /|
         /          / |
        +----------+  |
        |          |  |
        |          |  +
        2.         | /
        |          |/
        +--- 3. ---+

    The first dimension holds the the 2D Arrays. Either 1 or 2 dimensions,
    depends on the calculation method (1 = Gauss-Seidel, 2 = Jacobi)
    The second dimension holds all the rows (rows + interlines). We will
    iterate over those first. This value changes if multiple processes/threads
    are calculating together.
    The third dimension holds the columns and are always N+1 long.
    */
    
    uint64_t i, j; /*  local variables for loops   */

    /* local variables, just copies of the arguments struct */
    uint64_t const N = arguments->N;
    uint64_t const N_rank = arguments->N_rank;

    arguments->M = allocateMemory(arguments->num_matrices * (N_rank) * (N + 1) * sizeof (double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof (double**));

    for (i = 0; i < arguments->num_matrices; i++) {
        arguments->Matrix[i] = allocateMemory(N_rank * sizeof (double*));

        for (j = 0; j < N_rank; j++) {
            arguments->Matrix[i][j] = arguments->M + (i * N_rank * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices(struct calculation_arguments* arguments, struct options const* options) {
    uint64_t g, i, j; /*  local variables for loops   */

    /* local variables, just copies of the arguments struct */
    uint64_t const N = arguments->N;
    uint64_t const N_rank = arguments->N_rank;
    double const h = arguments->h;
    double*** Matrix = arguments->Matrix;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++) {
        for (i = 0; i < N_rank; i++) {
            for (j = 0; j <= N; j++) {
                Matrix[g][i][j] = 0.0;
            }
        }
    }

    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0) {
        for (g = 0; g < arguments->num_matrices; g++) {
            /* interation walks only over the rank specific amount of lines/rows (N_rank) */
            /* to get the actual line an offset of row_start has to be added */
            for (i = 0; i < N_rank; i++) {
                if (options->rank == 0) {
                    /* the (real) first line get initalized from 1.0 ... 0.0 */
                    uint64_t j = 0;
                    for (j = 0; j <= N; j++) {
                        Matrix[g][0][j] = 1.0 - (h * j);
                    }
                }
                //kein else-if, da sonst bei einem Aufruf mit nur einem Prozess
                //die unterste Zeile nicht initialisiert wird
                if (options->rank == options->size - 1) {
                    /* the (real) last line get initalized from 0.0 ... 1.0 */
                    uint64_t j = 0;
                    for (j = 0; j <= N; j++) {
                        Matrix[g][N_rank - 1][j] = h * j;
                    }
                }

                /* initalize all right and left boarder varaibles */
                /* start_row has to be taken in the calculation */
                Matrix[g][i][0] = 1.0 - (h * (i + arguments->row_start));
                Matrix[g][i][N] = h * (i + arguments->row_start);
            }

            /* the Matrix fields should have already the values,    */
            /* just to make sure, that there are no rounding errors */
            Matrix[g][N_rank - 1][0] = 0.0;
            Matrix[g][0][N] = 0.0;
        }
    }
}



/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics(struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options) {
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof (double) * arguments->num_matrices / 1024.0 / 1024.0);
    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL) {
        printf("Gauss-Seidel");
    } else if (options->method == METH_JACOBI) {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %" PRIu64 "\n", options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0) {
        printf("f(x,y) = 0");
    } else if (options->inf_func == FUNC_FPISIN) {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Terminierung:       ");

    if (options->termination == TERM_PREC) {
        printf("Hinreichende Genaugkeit");
    } else if (options->termination == TERM_ITER) {
        printf("Anzahl der Iterationen");
    }

    printf("\n");
    printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
    printf("Norm des Fehlers:   %e\n", results->stat_precision);
    printf("\n");
}


/* ************************************************************************ */
/* calculate: solves the equation for Gauss-Seidel                          */
/* ************************************************************************ */
static
void
MPI_calculateGaussSeidel(struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options) {
    //MPI_Request up[2], down[2]; /* local varible to coordinate MPI requests */
	
	int i, j; /* local variables for loops  */
    int m1, m2; /* used as indices for old and new matrices       */
    double star; /* four times center value minus 4 neigh.b values */
    double residuum; /* residuum of current iteration                  */
    double maxresiduum; /* maximum residuum value of a slave in iteration */

     /* get local copy of variable form the structure */
    int const N = arguments->N; 
    int const N_rank = arguments->N_rank;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;
	int iteration_first = 1; //TRUE, if this is the first iteration

    /* initialize m1 and m2 */
    m1 = 0;
    m2 = 0;

    if (options->inf_func == FUNC_FPISIN) {
        pih = PI * h;
        fpisin = TWO_PI_SQUARE * h * h;
    }
	
	INFO_STR(options->rank, options->size, "MPI_calculateGaussSeidel START");

    while (term_iteration > 0) {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In = arguments->Matrix[m2];

        maxresiduum = 0;
		
		//Check first for termination signals (in case of termination after precision)
		if (options->termination == TERM_PREC) {
			
			if (!iteration_first && options->rank < options->size-1) {
				//Have we received a abort signal?
				
				//MPI_Iprobe(MPI_ANY_SOURCE, MPI_TAG_ABORT, MPI_COMM_WORLD, (void*)&residuum, NULL);
				INFO_STR(options->rank, options->size, "RECV Check for termination signal from ABOVE v");
                MPI_Recv(&term_iteration, 1, MPI_INT, MPI_ANY_SOURCE, MPI_TAG_ABORT, MPI_COMM_WORLD, NULL);
				
				if(term_iteration == 0) {
					//Yes, we have received one!
					INFO_STR(options->rank, options->size, "Preparing termination");
				}
				
				//Inform the next rank about it!
				if (options->rank < options->size-1 -1) {				
					//Send the termination signal
					INFO_STR(options->rank, options->size, "SEND termination signal to next rank v");
					MPI_Send(&term_iteration, 1, MPI_INT, options->rank+1, MPI_TAG_ABORT, MPI_COMM_WORLD);
				}
					
				if(term_iteration == 0) {
					//Yes, we have received one!
					//Stop the loop!
					break;
				}
				
				INFO_STR(options->rank, options->size, "NO termination signal from ABOVE v");
			
			
				//receive maxresiduum(s) of previous rank(s)
				if(0 < options->rank) {
					INFO_STR(options->rank, options->size, "RECV maxresiduum from ABOVE v");
					MPI_Recv(&maxresiduum, 1, MPI_INT, options->rank-1, MPI_TAG_PRECISION, MPI_COMM_WORLD, NULL);
				}
			}
		}
		
		if(0 < options->rank) {
			//Get the line from above.
			INFO_STR(options->rank, options->size, "RECV line from ABOVE v");
			MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, options->rank - 1, MPI_TAG_DOWN, MPI_COMM_WORLD, NULL);
			//MPI_Irecv(Matrix_Out[0], N + 1, MPI_DOUBLE, options->rank - 1, MPI_TAG_DOWN, MPI_COMM_WORLD, &down[1]);
			INFO_STR(options->rank, options->size, "WAIT RECV line from ABOVE v");
            //MPI_Wait(&down[1], NULL);
			INFO_STR(options->rank, options->size, "WAIT DONE line from ABOVE v");
            INFO_DATA(options->rank, options->size, Matrix_Out[0] , N + 1);
		}
		if (!iteration_first && options->rank < options->size-1) {
			//Get the (last) line from below.
			INFO_STR(options->rank, options->size, "RECV line from BELOW ^");
			//MPI_Irecv(Matrix_Out[N_rank - 1], N + 1, MPI_DOUBLE, options->rank + 1, MPI_TAG_UP, MPI_COMM_WORLD, &up[1]);
			MPI_Recv(Matrix_Out[N_rank - 1], N + 1, MPI_DOUBLE, options->rank + 1, MPI_TAG_UP, MPI_COMM_WORLD, NULL);
		}
		
		INFO_STR(options->rank, options->size, "Calculate START");
		
        /* over all rows */
        for (i = 1; i < N_rank-1; i++) {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN) {
                /* add offset */
                fpisin_i = fpisin * sin((double) (i + arguments->row_start) * pih);
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                /* Matrix access is still in the i-index, no offset */
                star = Matrix_In[i][j] - (0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]));
                
				residuum = -star;
                if (options->inf_func == FUNC_FPISIN) {
					residuum = (fpisin_i * sin((double)(j) * pih)) - star;
                }

                Matrix_Out[i][j] = Matrix_In[i][j] + residuum;
				
				residuum = (residuum < 0) ? -residuum : residuum;
				maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
            }
			
			//Special cases, to improve speed (and make gauss-seidel work)			
			if(i==1) {
				//First line was calculated.
				if(0 < options->rank) {
					//Now send fast to the previous rank (last line of the previous rank)
					INFO_STR(options->rank, options->size, "SEND line ABOVE ^");
					//MPI_Isend(Matrix_Out[1], N + 1, MPI_DOUBLE, options->rank - 1, MPI_TAG_UP, MPI_COMM_WORLD, &up[0]);
					MPI_Send(Matrix_Out[1], N + 1, MPI_DOUBLE, options->rank - 1, MPI_TAG_UP, MPI_COMM_WORLD);
					INFO_DATA(options->rank, options->size, Matrix_Out[1] , N + 1);
				}
			} else if(i == N_rank-2) {
				//Before the last line gets calculated
				if (!iteration_first && options->rank < options->size-1) {
					//Ensure that the line from the rank below arrived.
					INFO_STR(options->rank, options->size, "WAIT RECV line from BELOW ^");
					//MPI_Wait(&up[1], NULL);
					INFO_STR(options->rank, options->size, "WAIT DONE line from BELOW ^");
                    INFO_DATA(options->rank, options->size, Matrix_Out[N_rank - 1] , N + 1);
				}
			}
        }
		
		INFO_STR(options->rank, options->size, "Calculate END");
		
		/* Communicate the line to the other ranks */
        if (options->rank < options->size-1) {
			if(!iteration_first) {
				INFO_STR(options->rank, options->size, "WAIT SEND line BELOW v DONE");
				//MPI_Wait(&down[0], NULL);
				INFO_STR(options->rank, options->size, "WAIT DONE line BELOW v DONE");
			}
			INFO_STR(options->rank, options->size, "SEND line BELOW v");
            //Send the last line to next rank
			//MPI_Isend(Matrix_Out[N_rank - 2], N + 1, MPI_DOUBLE, options->rank + 1, MPI_TAG_DOWN, MPI_COMM_WORLD, &down[0]);
			MPI_Send(Matrix_Out[N_rank - 2], N + 1, MPI_DOUBLE, options->rank + 1, MPI_TAG_DOWN, MPI_COMM_WORLD);
			INFO_DATA(options->rank, options->size, Matrix_Out[N_rank - 2] , N + 1);
        }
        if(0 < options->rank) {
			//Ensure, that the data line was sent
			INFO_STR(options->rank, options->size, "WAIT SEND line ABOVE ^");
            //MPI_Wait(&up[0], NULL);
			INFO_STR(options->rank, options->size, "WAIT DONE line ABOVE ^");
        }
		
        results->stat_iteration++;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation, depending on termination method */
        if (options->termination == TERM_PREC) {
			if (options->rank == options->size-1) {
				INFO_DOUBLESTR(options->rank, options->size, "maxresiduum: %f", maxresiduum);
				if (maxresiduum < options->term_precision) {
					//Send the abort signal to first rank.
                    term_iteration = 0;
					INFO_STR(options->rank, options->size, "SEND termination signal to next rank v (initial)");
					MPI_Send(&term_iteration, 1, MPI_INT, 0, MPI_TAG_ABORT, MPI_COMM_WORLD);
				} else {
					//Send an info (no termination) to the first rank
					INFO_STR(options->rank, options->size, "SEND NO termination signal to next rank v (initial)");
                    MPI_Send(&term_iteration, 1, MPI_INT, 0, MPI_TAG_ABORT, MPI_COMM_WORLD);
                }
			}
			
			//send maxresiduum (so far) to next rank
			if (options->rank < options->size-1) {
				INFO_STR(options->rank, options->size, "SEND maxresiduum to next rank v");
				MPI_Send(&maxresiduum, 1, MPI_INT, options->rank+1, MPI_TAG_PRECISION, MPI_COMM_WORLD);
			}
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
        }
		
		iteration_first = 0;
    }
	
	INFO_STR(options->rank, options->size, "REDUCE");
	MPI_Allreduce(&maxresiduum, &(results->stat_precision), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	INFO_STR(options->rank, options->size, "REDUCE DONE");

    results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation for Jacobi                                */
/* ************************************************************************ */
static
void
MPI_calculateJacobi(struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options) {
    MPI_Request up[2], down[2]; /* local varible to coordinate MPI requests */

    int i, j; /* local variables for loops  */
    int m1, m2; /* used as indices for old and new matrices       */
    double star; /* four times center value minus 4 neigh.b values */
    double residuum; /* residuum of current iteration                  */
    double maxresiduum; /* maximum residuum value of a slave in iteration */

    /* get local copy of variable form the structure */
    int const N = arguments->N; 
    int const N_rank = arguments->N_rank;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 for Jacobi */
    m1 = 0;
    m2 = 1;

    if (options->inf_func == FUNC_FPISIN) {
        pih = PI * h;
        fpisin = TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0) {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In = arguments->Matrix[m2];

        maxresiduum = 0;

        /* over all rows */
        /* interation walks only over the rank specific amount of lines/rows (N_rank) */
        /* to get the actual line an offset of row_start has to be added */
        for (i = 1; i < N_rank - 1; i++) {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN) {
                /* add offset */
                fpisin_i = fpisin * sin(pih * (double) (i + arguments->row_start));
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                /* Matrix access is still in the i-index, no offset */
				
                star = Matrix_In[i][j] - (0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]));
                
				residuum = -star;
                if (options->inf_func == FUNC_FPISIN) {
					residuum = (fpisin_i * sin((double)(j) * pih)) - star;
                }

                Matrix_Out[i][j] = Matrix_In[i][j] + residuum;
				
				residuum = (residuum < 0) ? -residuum : residuum;
				maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
            }
        }

        /* Communicate the line to the other ranks */
        if (options->rank != options->size - 1) {
            /* Last line has to be send to the rank below */
            MPI_Isend(Matrix_Out[N_rank - 2], N + 1, MPI_DOUBLE, options->rank + 1, 0, MPI_COMM_WORLD, &down[0]);
            MPI_Irecv(Matrix_Out[N_rank - 1], N + 1, MPI_DOUBLE, options->rank + 1, 0, MPI_COMM_WORLD, &up[1]);
        }
        if (options->rank != 0) {
            /* the first line has to be send to the rank above */
            MPI_Isend(Matrix_Out[1], N + 1, MPI_DOUBLE, options->rank - 1, 0, MPI_COMM_WORLD, &up[0]);
            MPI_Irecv(Matrix_Out[0], N + 1, MPI_DOUBLE, options->rank - 1, 0, MPI_COMM_WORLD, &down[1]);
        }

        /* wait for the request to finish */
        if (options->rank != 0) {
            MPI_Wait(&down[1], NULL);
            MPI_Wait(&up[0], NULL);
        }
        if (options->rank != options->size - 1) {
            MPI_Wait(&down[0], NULL);
            MPI_Wait(&up[1], NULL);
        }

        /* Build the maxresiduum over all ranks */
        MPI_Allreduce(&maxresiduum, &(results->stat_precision), 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        results->stat_iteration++;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation, depending on termination method */
        if (options->termination == TERM_PREC) {
            if (results->stat_precision < options->term_precision) {
                term_iteration = 0;
            }
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    results->m = m2;
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
DisplayMatrix(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options) {
    int x, y;

    double** Matrix = arguments->Matrix[results->m];

    int const interlines = options->interlines;

    printf("Matrix:\n");

    for (y = 0; y < 9; y++) {
        for (x = 0; x < 9; x++) {
            printf("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
        }

        printf("\n");
    }

    fflush(stdout);
}

/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
MPI_DisplayMatrix(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to) {
    int const elements = 8 * options->interlines + 9;

    int x, y;
    double** Matrix = arguments->Matrix[results->m];
    MPI_Status status;

    /* first line belongs to rank 0 */
    if (rank == 0)
        from--;

    /* last line belongs to rank size - 1 */
    if (rank + 1 == size)
        to++;

    if (rank == 0)
        printf("Matrix:\n");

    for (y = 0; y < 9; y++) {
        int line = y * (options->interlines + 1);

        if (rank == 0) {
            /* check whether this line belongs to rank 0 */
            if (line < from || line > to) {
                /* use the tag to receive the lines in the correct order
                 * the line is stored in Matrix[0], because we do not need it anymore */
                MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
            }
        } else {
            if (line >= from && line <= to) {
                /* if the line belongs to this process, send it to rank 0
                 * (line - from + 1) is used to calculate the correct local address */
                MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
            }
        }

        if (rank == 0) {
            for (x = 0; x < 9; x++) {
                int col = x * (options->interlines + 1);

                if (line >= from && line <= to) {
                    /* this line belongs to rank 0 */
                    printf("%7.4f", Matrix[line][col]);
                } else {
                    /* this line belongs to another rank and was received above */
                    printf("%7.4f", Matrix[0][col]);
                }
            }

            printf("\n");
        }
    }

    fflush(stdout);
}


/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main(int argc, char** argv) {
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    //Preinitalize MPI related variables
    options.rank = -1;
    options.size = -1;

    //Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &(options.rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(options.size));

    /* get parameters */
    AskParams(&options, argc, argv);

    /*  get and initialize variables and matrices  */
    initVariables(&arguments, &results, &options);
    allocateMatrices(&arguments);
    initMatrices(&arguments, &options);

    gettimeofday(&start_time, NULL); /*  start timer  */
    if (options.method == METH_JACOBI) {
        /*  solve the equation for Jacobi  */
        MPI_calculateJacobi(&arguments, &results, &options); 
    } else {
        /*  solve the equation for Gauss-Seidel  */
        MPI_calculateGaussSeidel(&arguments, &results, &options); 
    }
    gettimeofday(&comp_time, NULL); /*  stop timer  */

    
    /* print the statistics (once) */
    if (options.rank <= 0) {
        displayStatistics(&arguments, &results, &options);
    }

    /* print the results */
    MPI_DisplayMatrix(&arguments, &results, &options, options.rank, options.size, arguments.row_start + 1, arguments.row_end - 1);

    //Clean up MPI
    MPI_Finalize();

    //Clean up program
    freeMatrices(&arguments); /*  free memory  */

    return 0;
}
