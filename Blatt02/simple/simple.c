/*
** simple err-demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *
mistake1 ()
{
  /*Lokale Variablen können nicht zurückgegeben werden. Somit muss
    Speicher allokiert werden (malloc)*/
  int *buf = malloc(sizeof(int) * 6);
  
  /*Array initializieren*/
  buf[0] = 1;
  int i;
  for(i=1; i<6; i++) {
      buf[i] = i;
  }
  
  return buf;
}

int *
mistake2 ()
{
  /*Für ein int-Array müssen auch Speicherblöcke in der Größe von int's
    allokiert werden. Nur sizeof(char) reicht nicht*/
  int *buf = malloc (sizeof (int) * 4);
  buf[1] = 2;
  return buf;
}

int *
mistake3 ()
{
  /*Auch über einen Pointertrick können keine Variablen bzw. Addressen
    von lokalen Variablen rausgegeben werden.*/
    
  /*int mistake2_ = 0;
  int *buf = (int *) &mistake2;*/
    
  int *buf = malloc (sizeof (int) * 4);  
  buf[0] = 3;
  return buf;
}

int *
mistake4 ()
{
  int *buf = malloc (sizeof (char) * 4); /* 4*char == 1*int */
  buf[0] = 4;
  /*Der Speicher darf nicht frühzeitig freigegeben werden.
    Nach Freigabe kann nicht mehr auf den Wert an der Addresse zugegriffen werden.
    free (buf);*/
  
  return buf;
}

int
main (void)
{   
  /* Modifizieren Sie die folgenden Zeilen nicht */
  int *p[4] = { &mistake1 ()[1], &mistake2 ()[1], mistake3 (), mistake4 () };

  printf ("1 %d\n", *p[0]);
  printf ("2 %d\n", *p[1]);
  printf ("3 %d\n", *p[2]);
  printf ("4 %d\n", *p[3]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free (p[0]-1);
  free (p[1]-1);
  free (p[2]);
  free (p[3]);

  return 0;
}
