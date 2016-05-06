AR=ar
RANLIB=ranlib
CC=mpicc
CFLAGS=-O2 -I./ -std=gnu99 -I$(HOME)/loc-install/include  -fopenmp -I$(HOME)/QuaHoG/lib -I$(HOME)/MG4QCD/include
LDFLAGS=-L$(HOME)/loc-install/lib -L$(HOME)/MG4QCD/lib -L$(HOME)/QuaHoG/lib
LIBS=-llime -fopenmp -lqhg -lmg4qcd -lhdf5 -lmxml

all: main

QHG=$(HOME)/QuaHoG/lib/libqhg.a
MG4QCD=$(HOME)/MG4QCD/lib/libmg4qcd.a

# A non-exaustive list of dependencies
main: main.o parser.o mg4qcd_interface.o $(QHG) $(MG4QCD) 
	$(CC) $(CFLAGS) -o main main.o parser.o mg4qcd_interface.o $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o
