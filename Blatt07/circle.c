#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int*
init (int N)
{
  //todo
  int* buf = malloc(sizeof(int) * N);

  srand(time(NULL));

  for (int i = 0; i < N; i++)
  {
    buf[i] = rand() % 25; //do not modify %25
  }

  return buf;
}

int*
circle (int* buf)
{
  //todo
  return buf;
}

int
main (int argc, char** argv)
{
  char arg[256];
  int N;
  int rank;
  int* buf;

  if (argc < 2)
  {
    printf("Arguments error\n");
    return EXIT_FAILURE;
  }

  sscanf(argv[1], "%s", arg);

  //array length
  N = atoi(arg);
  buf = init(N);

  //todo myrank
  rank = 0;

  printf("\nBEFORE\n");

  for (int i = 0; i < N; i++)
  {
    printf ("rank %d: %d\n", rank, buf[i]);
  }

  circle(buf);

  printf("\nAFTER\n");

  for (int j = 0; j < N; j++)
  {
    printf ("rank %d: %d\n", rank, buf[j]);
  }

  return EXIT_SUCCESS;
}
