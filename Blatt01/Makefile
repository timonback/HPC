all: map pointer
	
map: map.o
	gcc -o map map.o
	
map.o: map.c
	gcc -c map.c
	
pointer: pointer.o
	gcc -o pointer pointer.o
	
pointer.o: pointer.c
	gcc -c pointer.c

clean:
	rm -f map
	rm -f pointer
	rm -f *.o
	rm -f *~