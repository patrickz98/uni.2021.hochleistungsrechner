# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff partdiff-posix
OBJS = partdiff.o askparams.o

# Targets ...
all: partdiff partdiff-posix

partdiff-posix: $(OBJS) Makefile
	$(CC) $(LFLAGS) -pthread -DPOSIX -o $@ $(OBJS) $(LIBS)

partdiff: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff.o: partdiff.c Makefile
	$(CC) -c $(CFLAGS) partdiff.c -o partdiff.o

askparams.o: askparams.c Makefile
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
