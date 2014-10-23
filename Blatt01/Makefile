all: pointer map

CC: gcc
CFLAGS: -Wall -Wextra

#define objects that are made from .c files
OBJ: pointer.o map.o

#make clean a phony target so that it will always be created
.PHONY: clean

#make target
%: %.o
	$(CC) $(CFLAGS) -o $@ $<

#make objects
%.o: %.c
	$(CC) $(CFLAGS) -c $<

#clean up
clean:
	rm -f pointer
	rm -f map
	rm -f *.o
	rm -f *~
