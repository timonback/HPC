#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <mpi.h>

#define UPPER_PROCESS(currentRank, rankSize) (currentRank+1 < rankSize ? currentRank+1 : 0)
#define LOWER_PROCESS(currentRank, rankSize) (currentRank   > 0 ? currentRank-1 : rankSize-1)

#define UNUSED(arg) (void)arg;

int*
init (int N, int rank)
{
    int* buf = (int*)malloc(sizeof(int) * N);

    srand(time(NULL) + (15*rank) );

    for (int i = 0; i < N; i++)
    {
        buf[i] = rand() % 25; //do not modify %25
    }

    return buf;
}

int*
circle (int* buf, int rank, int size, int* data_my_width, int* data_rank_width)
{
    MPI_Send(buf, *data_my_width, MPI_INT, UPPER_PROCESS(rank, size), 0, MPI_COMM_WORLD);
    MPI_Recv(buf, *data_rank_width, MPI_INT, LOWER_PROCESS(rank, size), 0, MPI_COMM_WORLD, NULL);
	
	(*data_my_width) = *data_rank_width;
	
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
    int data_my_start = data_width * rank;
    int data_my_end = data_width * (rank + 1);
    int data_my_width = data_my_end - data_my_start;
	int data_previous_rank_width = 0;

    buf = init(data_my_width+1, rank);
	
	//Save the termination value (just interesting for the last process)
    if(rank == 0) {
		MPI_Send(buf, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD);
	} else if (rank == size - 1) {
		MPI_Recv(&term_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
	}
	
	//START UP OUTPUT
    if (rank == 0)
    {
        printf("\nBEFORE\n");

        for (int i = 0; i < data_my_width; i++)
        {
           printf ("rank %d: %d\n", rank, buf[i]);
        }
		
		int buf_temp[data_my_width+1];
        for (int i = 1; i < size; i++)
        {
            int rank_realWidth;
            MPI_Recv(&rank_realWidth, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            MPI_Recv(buf_temp, rank_realWidth, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            
            for (int j = 0; j < rank_realWidth; j++)
            {
                printf ("rank %d: %d\n", i, buf_temp[j]);
            }
        }
    }
    else
    {
        MPI_Send(&data_my_width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(buf, data_my_width, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
	//Synchronize all processes.
	//So the processes will continue, AFTER the main process has printed everything.
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}

	
	
	//Synchronize all processes.
	//So the processes will continue, AFTER the main process has informed everyone about the termination value
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}

	
	//Circle loop
	int circle_run_loop = 1;
    while(circle_run_loop)
    {
        MPI_Send(&data_my_width, 1, MPI_INT, UPPER_PROCESS(rank, size), 0, MPI_COMM_WORLD);
        MPI_Recv(&data_previous_rank_width, 1, MPI_INT, LOWER_PROCESS(rank, size), 0, MPI_COMM_WORLD, NULL);

        buf = circle(buf, rank, size, &data_my_width, &data_previous_rank_width);
		
        if (rank == size - 1)
        {
            circle_run_loop = (buf[0] != term_value);
            for (int i = 0; i < size - 1; i++)
            {
                MPI_Send(&circle_run_loop, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Recv(&circle_run_loop, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, NULL);
		}
		
		//Synchronize all processes.
		//So the processes will continue, AFTER all process know, if the termination will/has been reached
		if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
			printf("Error in MPI_Barrier\n");
		}
		
    }
	
	
	
	//END OUTPUT
    if (rank == 0)
    {
        printf("\nAFTER\n");

        for (int i = 0; i < data_my_width; i++)
        {
            printf ("rank %d: %d\n", rank, buf[i]);
        }
		
		int buf_temp[data_my_width+1];
        for (int i = 1; i < size; i++)
        {
            int rank_realWidth;
            MPI_Recv(&rank_realWidth, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            MPI_Recv(buf_temp, rank_realWidth, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
			
            for (int j = 0; j < rank_realWidth; j++)
            {
                printf ("rank %d: %d\n", i, buf_temp[j]);
            }
        }
    }
    else
    {
        MPI_Send(&data_my_width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(buf, data_my_width, MPI_INT, 0, 0, MPI_COMM_WORLD);
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
