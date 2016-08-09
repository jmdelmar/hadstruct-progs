AR=ar
RANLIB=ranlib
CC=mpicc
CFLAGS=-O2 -I./ -std=gnu99 -I$(HOME)/loc-install/include  -fopenmp -I$(HOME)/QuaHoG/lib -I$(HOME)/DDalphaAMG/include
LDFLAGS=-L$(HOME)/loc-install/lib -L$(HOME)/DDalphaAMG/lib -L$(HOME)/QuaHoG/lib
LIBS=-llime -fopenmp -lqhg -lDDalphaAMG -lhdf5 -lmxml
all: main

QHG=$(HOME)/QuaHoG/lib/libqhg.a
MG=$(HOME)/MG4QCD/lib/libDDalphaAMG.a

# A non-exaustive list of dependencies
main: main.o parser.o mg_interface.o $(QHG) $(MG4QCD) 
	$(CC) $(CFLAGS) -o main main.o parser.o mg_interface.o $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o
