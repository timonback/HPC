#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mpi.h>

#define UPPER_PROCESS(currentRank, rankSize) (currentRank+1 < rankSize ? currentRank+1 : 0)
#define LOWER_PROCESS(currentRank, rankSize) (currentRank   > 0 ? currentRank-1 : rankSize-1)

#define UNUSED(arg) (void)arg;

int*
init (int N, int rank)
{
    int* buf = malloc(sizeof(int) * N);

    srand(time(NULL) + (15*rank) );

    for (int i = 0; i < N; i++)
    {
        buf[i] = rand() % 25; //do not modify %25
    }

    return buf;
}

int*
circle (int* buf, int step_width, int previous_width, int rank, int size)
{
    MPI_Send(buf, step_width, MPI_INT, UPPER_PROCESS(rank, size), 0, MPI_COMM_WORLD);

    buf = realloc(buf, previous_width * sizeof(int));
    
    MPI_Recv(buf, previous_width, MPI_INT, LOWER_PROCESS(rank, size), 0, MPI_COMM_WORLD, NULL);
	
	//Synchronize all processes.
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}

    return buf;
}

void usage(int argc, char* argv[])
{
	UNUSED(argc);
	
	printf("This program generates an N many random integers array and rotates them in a circle. Every possible core of the system will be used.\n");
	printf("It will loop as long as the first item of the first process at the initialization and first item of last process at that specific iteration will NOT match. ");
	printf("So, maximum count of iterations will be NUM_CORES-1");
	printf("Usage: %s N\n", argv[0]);
	printf("Example: %s 10\n", argv[0]);
}

int
main (int argc, char** argv)
{
    char arg[256];
    int N;
    int* buf;
    int rank, size;
    int term_value;
    int compare_value;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2)
    {
        if (rank == 0)
        {
            printf("Arguments error. Please specify all required arguments\n");
			usage(argc, argv);
        }
    
        MPI_Finalize();
        
        return EXIT_FAILURE;
    }

    sscanf(argv[1], "%s", arg);

    //array length
    N = atoi(arg);
	
	if(N < size) {
		if (rank == 0)
        {
            printf("Your system has %i cores. The data, which will be circled, must be larger than the number of cores (data>%i)\n", size, size);
			usage(argc, argv);
        }
    
        MPI_Finalize();
        
        return EXIT_FAILURE;
	}

    double data_width = (double)N / (double)size;
    int data_start = data_width * rank;
    int data_end = data_width * (rank + 1);
    int data_realWidth = data_end - data_start;

    buf = init(data_realWidth, rank);

	
	//START UP OUTPUT
    if (rank == 0)
    {
        printf("\nBEFORE\n");

        for (int i = 0; i < data_realWidth; i++)
        {
           printf ("rank %d: %d\n", rank, buf[i]);
        }

        for (int i = 1; i < size; i++)
        {
            int rank_realWidth;
            MPI_Recv(&rank_realWidth, 1, MPI_INT, i, 5, MPI_COMM_WORLD, NULL);
            int buf_temp[rank_realWidth];

            MPI_Recv(buf_temp, rank_realWidth, MPI_INT, i, 1, MPI_COMM_WORLD, NULL);
            
            for (int j = 0; j < rank_realWidth; j++)
            {
                printf ("rank %d: %d\n", i, buf_temp[j]);
            }
        }
    }
    else
    {
        MPI_Send(&data_realWidth, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
        MPI_Send(buf, data_realWidth, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }
	//Synchronize all processes.
	//So the processes will continue, AFTER the main process has printed everything.
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}

	
	//Send (start) termination value to all
    if (rank == 0)
    {
        term_value = buf[0];
        for (int i = 1; i < size; i++)
        {
            MPI_Send(&term_value, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&term_value, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, NULL);
    }
	//Synchronize all processes.
	//So the processes will continue, AFTER the main process has informed everyone about the termination value
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}

	
	//Circle loop
    int previous_width;
    do
    {
        MPI_Send(&data_realWidth, 1, MPI_INT, UPPER_PROCESS(rank, size), 7, MPI_COMM_WORLD);
        MPI_Recv(&previous_width, 1, MPI_INT, LOWER_PROCESS(rank, size), 7, MPI_COMM_WORLD, NULL);

        circle(buf, data_realWidth, previous_width, rank, size);
		
        if (rank == size - 1)
        {
            compare_value = buf[0];
            for (int i = 0; i < size - 1; i++)
            {
                MPI_Send(&compare_value, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&compare_value, 1, MPI_INT, size - 1, 3, MPI_COMM_WORLD, NULL);
		}
		
		//Synchronize all processes.
		//So the processes will continue, AFTER all process know, if the termination will/has been reached
		if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
			printf("Error in MPI_Barrier\n");
		}
		
    } while (term_value != compare_value);
	
	
	
	//END OUTPUT
    if (rank == 0)
    {
        printf("\nAFTER\n");

        for (int i = 0; i < data_realWidth; i++)
        {
            printf ("rank %d: %d\n", rank, buf[i]);
        }

        for (int i = 1; i < size; i++)
        {
            int rank_realWidth;
            MPI_Recv(&rank_realWidth, 1, MPI_INT, i, 6, MPI_COMM_WORLD, NULL);
            int buf_temp[rank_realWidth];

            MPI_Recv(buf_temp, rank_realWidth, MPI_INT, i, 4, MPI_COMM_WORLD, NULL);
            for (int j = 0; j < rank_realWidth; j++)
            {
                printf ("rank %d: %d\n", i, buf_temp[j]);
            }
        }
    }
    else
    {
        MPI_Send(&data_realWidth, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
        MPI_Send(buf, data_realWidth, MPI_INT, 0, 4, MPI_COMM_WORLD);
    }
	//Synchronize all processes.
	//So the processes will continue, AFTER the main process has printed everything.
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}

	//Cleanup
    free(buf);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
