AR=ar
RANLIB=ranlib
CC=cc
CFLAGS=-O2 -I./ -I$(HOME)/loc-install/include -I$(HOME)/QuaHoG/lib -I$(HOME)/MG4QCD/include
LDFLAGS=-L$(HOME)/loc-install/lib -L$(HOME)/MG4QCD/lib -L$(HOME)/QuaHoG/lib
LIBS=-lqhg -lmg4qcd -lhdf5 -lmxml -llime

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
