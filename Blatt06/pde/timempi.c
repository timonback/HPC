/*
echo $(hostname): $(($(date +%s%N)/1000000))

exit
*/

/*
	srun -N [Anzahl der Nodes] mpiexec \-n [Anzahl der Prozesse] ./[Programm]
*/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	int i=0;
	char* host = (char*)malloc(sizeof(char)*255);
	gethostname(host, 255);
	for(i=0; i<argc; i++) {
		printf("%d:%i: %s %s\n", 0, i, host, argv[i]);
	}

	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	printf("Hello World from process %d of %d\n", rank, size);
	MPI_Finalize();
	return 0;
}