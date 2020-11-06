all: map pointer

clean: map.clean pointer.clean

map: map.o
	gcc -o map map.o

map.o: map.c
	gcc -c map.c

map.clean:
	rm -f map
	rm -f *.o
	rm -f *~

pointer: pointer.o
	gcc -o pointer pointer.o
	
pointer.o: pointer.c
	gcc -c pointer.c

pointer.clean:
	rm -f pointer
	rm -f *.o
	rm -f *~
