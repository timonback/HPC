#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mpi.h>

int*
init (int N, int rank)
{
    int* buf = malloc(sizeof(int) * N);

    srand(time(NULL) + rank);

    for (int i = 0; i < N; i++)
    {
        buf[i] = rand() % 25; //do not modify %25
    }

    return buf;
}

int*
circle (int* buf, int step_width, int previous_width, int rank, int size)
{
    if (rank < size - 1)
    {
        MPI_Send(buf, step_width, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Send(buf, step_width, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    buf = realloc(buf, previous_width * sizeof(int));
    
    if (rank > 0)
    {
        MPI_Recv(buf, step_width + 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, NULL);
    }
    else
    {
        MPI_Recv(buf, step_width + 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, NULL);
    }

    return buf;
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
            printf("Arguments error\n");
            printf("This program generates an N many random integers array and rotates them in a circle. Every possible core of the system will be used.\n");
            printf("It will loop as long as the first item of the first process at the initialization and first item of last process at that specific iteration will NOT match. ");
            printf("So, maximum count of iterations will be NUM_CORES-1");
            printf("Usage: %s N\n", argv[0]);
            printf("Example: %s 10\n", argv[0]);
        }
    
        MPI_Finalize();
        
        return EXIT_FAILURE;
    }

    sscanf(argv[1], "%s", arg);

    //array length
    N = atoi(arg);

    double data_width = (double)N / (double)size;
    int data_start = data_width * rank;
    int data_end = data_width * (rank + 1);
    int data_realWidth = data_end - data_start;

    buf = init(data_realWidth, rank);

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

    int num_of_iterations = 0;

    do
    {
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

        int previous_width;

        if(rank > 0)
        {
            MPI_Send(&data_realWidth, 1, MPI_INT, rank - 1, 7, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&data_realWidth, 1, MPI_INT, size - 1, 7, MPI_COMM_WORLD);
        }

        if(rank < size - 1)
        {
            MPI_Recv(&previous_width, 1, MPI_INT, rank + 1, 7, MPI_COMM_WORLD, NULL);
        }
        else
        {
            MPI_Recv(&previous_width, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, NULL);
        }

        circle(buf, data_realWidth, previous_width, rank, size);

        num_of_iterations++;
    } while (term_value != compare_value && num_of_iterations < size - 1);

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

    free(buf);

    MPI_Finalize();

    return EXIT_SUCCESS;
}
