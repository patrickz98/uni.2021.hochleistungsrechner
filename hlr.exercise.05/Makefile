# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-posix
OBJS = partdiff.o askparams.o

# Targets ...
all: partdiff

partdiff: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $(TGTS) $(OBJS) $(LIBS)

partdiff.o: partdiff.c Makefile

partdiff-openmp.o: partdiff-openmp.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
	$(RM) partdiff-openmp

openmp: partdiff-openmp.o askparams.o
	$(CC) $(LFLAGS) -o partdiff-openmp partdiff-openmp.o askparams.o $(LIBS)

run: clean all
	# OMP_NUM_THREADS=12
	# ./partdiff-openmp 12 2 64 2 2 5120 > zzz.openmp.txt
	./partdiff-posix 12 2 64 2 2 5120 > zzz.posix.txt
	# ./partdiff-posix 12 2 64 2 2 100 > zzz.posix.100.txt
	# ./partdiff-posix 12 2 64 2 2 2 > zzz.posix.2.txt
