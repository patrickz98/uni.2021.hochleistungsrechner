CC = mpicc
CFLAGS = -pedantic -Wall -Wextra -Wno-unused-function -g

TGTS = timempi timempi2

all: clean $(TGTS)

timempi : timempi.o Makefile
	$(CC) $(CFLAGS) -o $@ timempi.o

timempi2 : timempi2.o Makefile
	$(CC) $(CFLAGS) -o $@ timempi2.o

timempi.o : timempi.c Makefile
	$(CC) -c $(CFLAGS) timempi.c -o timempi.o

timempi2.o : timempi2.c Makefile
	$(CC) -c $(CFLAGS) timempi2.c -o timempi2.o


clean:
	rm -f $(TGTS)
	rm -f *.o
	rm -f *~
