CC = mpicc
CFLAGS = -pedantic -Wall -Wextra -Wno-unused-function
LFLAGS = $(CFLAGS)

TGTS = circle

all: $(TGTS)

% : %.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f $(TGTS)
	rm -f *.o
	rm -f *~
