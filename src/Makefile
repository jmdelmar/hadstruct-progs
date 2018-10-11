AR=ar
RANLIB=ranlib
CC=mpiicc

#... Adjust according to compiler
CFLAGS=-O3 -D_GNU_SOURCE -I./ -std=gnu99 -fopenmp -xCORE-AVX2

#... Path to qhg.h
QHG_INC=-I$(HOME)/BDW/QuaHoG/lib
CFLAGS+=$(QHG_INC)

#... Path to lime.h
LIME_INC=-I$(HOME)/BDW/install/include
CFLAGS+=$(LIME_INC)

#... Path to mxml.h
MXML_INC=-I$(HOME)/BDW/install/include
CFLAGS+=$(MXML_INC)

#... Path to DDalphaAMG.h
MG_INC=-I$(HOME)/BDW/DDalphaAMG/include
CFLAGS+=$(MG_INC)

#... Path to hdf5.h
HDF5_INC=-I$(HDF5_BASE)/include
CFLAGS+=$(HDF5_INC)

#... Paths to libraries
LDFLAGS="-xCORE-AVX2"

#... Path to libqhg.a
QHG_LIB=-L$(HOME)/BDW/QuaHoG/lib
LDFLAGS+=$(QHG_LIB)

#... Path to liblime.a
LIME_LIB=-L$(HOME)/BDW/install/lib
LDFLAGS+=$(LIME_LIB)

#... Path to libmxml.a
MXML_LIB=-L$(HOME)/BDW/install/lib
LDFLAGS+=$(MXML_LIB)

#... Path to libDDalphaAMG.a
MG_LIB=-L$(HOME)/BDW/DDalphaAMG/lib
LDFLAGS+=$(MG_LIB)

#... Path to libhdf5.a
HDF5_LIB=-L$(HDF5_BASE)/lib
LDFLAGS+=$(HDF5_LIB)

#... Order may matter depending on compiler
LIBS=-fopenmp -lqhg -lDDalphaAMG -llime -lmxml -lhdf5

main: main.o mg_interface.o parser.o 
	$(CC) $(CFLAGS) -o $@ $< mg_interface.o parser.o $(LDFLAGS) $(LIBS)

mesons: main_mesons.o mg_interface.o parser.o 
	$(CC) $(CFLAGS) -o $@ $< mg_interface.o parser.o $(LDFLAGS) $(LIBS)


%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o
