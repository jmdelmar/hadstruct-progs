AR=ar
RANLIB=ranlib
CC=cc
CFLAGS=-D_GNU_SOURCE -I$(HOME)/tmLQCD/include -I$(HOME)/QuaHoG/lib
LDFLAGS=-L$(HOME)/tmLQCD/lib -L$(HOME)/loc-install/lib -L$(HOME)/QuaHoG/lib
LIBS=-lqhg -lwrapper -lhmc -lmonomial -loperator -lsolver -linit -llinalg -lhmc -lxchange -lrational -lio -lparpack -larpack -lm -llemon -llime 

all: main

# A non-exaustive list of dependencies
main: main.o $(HOME)/tmLQCD/lib/libwrapper.a $(HOME)/QuaHoG/lib/libqhg.a
	$(CC) $(CFLAGS) -o _main_ main.o $(LDFLAGS) $(LIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	$(RM) main.o
