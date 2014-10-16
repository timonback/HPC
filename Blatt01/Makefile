all: map
	
pointer: map.o
	gcc -o map map.o
	
pointer.o: map.c
	gcc -c map.c

clean:
	rm -f map
	rm -f *.o
	rm -f *~
