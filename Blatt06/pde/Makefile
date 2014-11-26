# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2
LFLAGS = $(CFLAGS)
LIBS   = -lmpi

OBJS_TIME = timempi.o
OBJS_TIME2 = timempi2.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: timempi timempi2

timempi: $(OBJS_TIME) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS_TIME) $(LIBS)
	
timempi: $(OBJS_TIME2) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS_TIME2) $(LIBS)

clean:
	$(RM) *.o *~

timempi.o: timempi.c Makefile
timempi2.o: timempi2.c Makefile