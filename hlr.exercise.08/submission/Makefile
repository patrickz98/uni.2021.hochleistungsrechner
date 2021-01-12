# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -pedantic -Wall -Wextra -Wno-unused-function -O3 -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff partdiff-par partdiff-par-hybrid
OBJS = partdiff.o askparams.o

# Targets ...
all: $(TGTS)

%: %.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ $*.o askparams.o $(LIBS)

partdiff.o:
	$(CC) -c $(CFLAGS) partdiff.c -o $*.o

partdiff-par.o: $(OBJS) Makefile
	$(CC) -c $(CFLAGS) partdiff.c -o $*.o -DPARALLEL

partdiff-par-hybrid.o: $(OBJS) Makefile
	$(CC) -c $(CFLAGS) partdiff.c -o $*.o -DPARALLEL -DHYBRID


askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c





clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
	$(RM) *.o
