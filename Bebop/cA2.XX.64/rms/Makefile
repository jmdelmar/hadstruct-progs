CC=mpicc
CFLAGS=-O2 -I./ -std=gnu99 -D_GNU_SOURCE -I$(HOME)/loc-install/include  -qopenmp -I$(HOME)/QuaHoG/lib
LDFLAGS=-L$(HOME)/loc-install/lib -L$(HOME)/QuaHoG/lib -L$(HDF5_BASE)/lib
LIBS=-qopenmp -lqhg -llime -lhdf5 -lmxml
all: main

QHG=$(HOME)/QuaHoG/lib/libqhg.a

# A non-exaustive list of dependencies
main: main.o $(QHG)
	$(CC) $(CFLAGS) -o main main.o  $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o
