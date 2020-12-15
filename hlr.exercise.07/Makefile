CC = mpicc
CFLAGS = -pedantic -Wall -Wextra -Wno-unused-function -g

TGTS = circle

all: clean $(TGTS)

circle : circle.o Makefile
	$(CC) $(CFLAGS) -o $@ circle.o

circle.o : circle.c Makefile
	$(CC) -c $(CFLAGS) circle.c -o circle.o

clean:
	rm -f $(TGTS)
	rm -f *.o
	rm -f *~

run: all
	mpirun -n 4 ./circle 17
