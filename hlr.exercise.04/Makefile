# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-openmp
OBJS = partdiff.o askparams.o

# Targets ...
all: partdiff

partdiff: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff.o: partdiff.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)

run: clean all
	# ./partdiff-openmp 1 2 64 2 2 5120 > zzz.out.txt
	# ./partdiff-openmp 1 2 512 2 2 5120 > zzz.out.txt
	# ./partdiff-openmp 1 2 64 2 2 10240 > zzz.out.64.10240.txt
	./partdiff-openmp 1 2 512 2 2 10240 > zzz.out.512.10240.txt
