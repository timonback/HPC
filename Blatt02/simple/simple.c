/*
** simple err-demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *
mistake1 ()
{
  //use no static
  static int buf[] = { 1, 1, 2, 3, 4, 5 };
  return buf;
}

int *
mistake2 ()
{
  int *buf = malloc (sizeof (int) * 4);
  buf[1] = 2;
  return buf;
}

int *
mistake3 ()
{
  //http://stackoverflow.com/questions/10038470/pointer-to-local-variable-in-c
  int mistake2_ = 0;
  int *buf = (int *) &(mistake2_);
  //buf-=3*sizeof(int);
  printf("%p\n", buf);
  printf("%p\n", &(buf[3]));
  printf("%p\n", &(buf[3])-3);
  *(&(buf[3])-3) = 3;
  return buf;
}

int *
mistake4 ()
{
  int *buf = malloc (sizeof (char) * 4); //4*char == 1*int
  buf[0] = 4;
  return buf;
}

int
main (void)
{
    /*
    int x = mistake1()[0];
    printf("value %d\n", x);
    printf("point; %d\n", *mistake1 ());
    printf("cvalue; %d\n", *(&mistake1 ()[0]));
    
    return 0;
    */
    
    /*ALIAS
     int p[4] = { mistake1 ()[1], mistake2 ()[1], *mistake3 (), *mistake4 () };

  printf ("1 %d\n", p[0]);
  printf ("2 %d\n", p[1]);
  printf ("3 %d\n", p[2]);
  printf ("4 %d\n", p[3]);
  
  return 0;
}
int remain() {*/
    
  /* Modifizieren Sie die folgenden Zeilen nicht */
  int *p[4] = { &mistake1 ()[1], &mistake2 ()[1], mistake3 (), mistake4 () };

  printf ("1 %d\n", *p[0]);
  printf ("2 %d\n", *p[1]);
  printf ("3 %d\n", *p[2]);
  printf ("4 %d\n", *p[3]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free (p[1]-1);			/* welcher Pointer war das doch gleich?, TODO: Fixme... :-) */
  free (p[3]);

  return 0;
}
