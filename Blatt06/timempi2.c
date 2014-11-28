#define _BSD_SOURCE //make gethostname workable without warning (non standard)

#include <stdio.h>  //printf
#include <string.h> //strcat, snprintf
#include <time.h>   //time, localtime, strftime

#include <unistd.h> //gethostname
#include <sys/time.h> //gettimeofday

#include <mpi.h>    //MPI_*

//constant for the buffer size
#define BUF_SIZE 255

int main(int argc, char* argv[])
{
	//initialize MPI
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	//buffer for receiving messages.
	char receive[BUF_SIZE];
	//integer variables to save the MicroSeconds.
	//usec, will be defaulted to the highest uSec possible. So that even unchanged uSecs
	//can be reduced with a usable output in lowest_usec.
	long int usec=999999, lowest_usec;
	
	if(rank == 0) {
		//Only the main process will receive data and print it.
		int i;
		for(i=1; i<size; ++i) {
			//Receive ...
			MPI_Recv(receive, BUF_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD, NULL);
			//... and print
			printf("%s\n", receive);
		}
	} else {
		//All other process will first compose a message string <Hostname: date.time.uSec>...
		char message[BUF_SIZE];
		
		//... containing the hostname ...
		if(gethostname(message, BUF_SIZE) != 0) {
			printf("Error in gethostname\n");
		}
		strcat(message, ": ");
		
		//... and the time ...
		// collect together the neccessary informations of date and time
		time_t date = time(NULL);
		struct tm *tmptr = localtime(&date);
		if(tmptr == NULL) {
			printf("Error in localtime\n");
		}
		struct timeval time;
		if(gettimeofday(&time, NULL) != 0) {
			printf("Error in gettimeofday\n");
		}
		// and print it in the message string
		strftime(message+strlen(message), BUF_SIZE-strlen(message), "%Y-%m-%d %H:%M:%S.", tmptr);
		// with the microseconds
		snprintf(message+strlen(message), BUF_SIZE-strlen(message), "%06ld", time.tv_usec);
		
		//... and finally send it.
		MPI_Send(message, BUF_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
		
		// save the microseconds of this process in usec.
		usec = time.tv_usec;
	}
	
	//Collect all the microseconds of the processes together and reduce them to the min value.
	MPI_Reduce(&usec, &lowest_usec, 1, MPI_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
	if(rank == 0) {
		//Only the main process will print the reduced microseconds value.
		printf("%06ld\n", lowest_usec);
	}
	
	//Synchronize all processes.
	//So the processes will continue, AFTER the main process has printed everything.
	if(MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
		printf("Error in MPI_Barrier\n");
	}
	
	//Processes may end now.
	printf("Rank %i beendet jetzt!\n", rank);
	MPI_Finalize();
	
	return 0;
}