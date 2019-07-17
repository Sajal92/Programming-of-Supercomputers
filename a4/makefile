CC = mpicc
SRC_GAUSS = gauss
SRC_MGEN = mgen

all: gauss mgen

gauss: $(SRC_GAUSS:%=%.o)
	$(CC) $(CFLAGS) -o gauss $(SRC_GAUSS:%=%.o)  -lm

mgen: $(SRC_MGEN:%=%.o)
	$(CC) $(CFLAGS) -o mgen $(SRC_MGEN:%=%.o)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	/bin/rm -f $(SRC_TEST:%=%.o) gauss mgen *.o
