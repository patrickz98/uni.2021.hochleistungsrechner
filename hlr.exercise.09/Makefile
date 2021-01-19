# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff
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
	mpirun -n 4 ./partdiff 1 1 64 1 2 8000
	# Referance f1
	# mpirun -n 4 ./partdiff 1 1 0 1 2 47
	# Referance f2
	# mpirun -n 4 ./partdiff 1 1 0 2 2 48
