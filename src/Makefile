AR=ar
RANLIB=ranlib
CC=mpicc

#... Adjust according to compiler
CFLAGS=-O2 -I./ -std=gnu99 -fopenmp -I$(HOME)/loc-install/include -I$(HOME)/DDalphaAMG/include

#... Path to qhg.h
QHG_INC=-I$(HOME)/QuaHoG/lib
CFLAGS+=$(QHG_INC)

#... Path to lime.h
LIME_INC=-I$(HOME)/loc-install/include
CFLAGS+=$(LIME_INC)

#... Path to mxml.h
MXML_INC=-I$(HOME)/loc-install/include
CFLAGS+=$(MXML_INC)

#... Path to DDalphaAMG.h
MG_INC=-I$(HOME)/DDalphaAMG/include
CFLAGS+=$(MG_INC)

#... Path to hdf5.h
HDF5_INC=
CFLAGS+=$(HDF5_INC)

#... Paths to libraries
LDFLAGS=

#... Path to qhg.h
QHG_LIB=-L$(HOME)/QuaHoG/lib
LDFLAGS+=$(QHG_LIB)

#... Path to lime.h
LIME_LIB=-L$(HOME)/loc-install/lib
LDFLAGS+=$(LIME_LIB)

#... Path to mxml.h
MXML_LIB=-L$(HOME)/loc-install/lib
LDFLAGS+=$(MXML_LIB)

#... Path to DDalphaAMG.h
MG_LIB=-L$(HOME)/DDalphaAMG/lib
LDFLAGS+=$(MG_LIB)

#... Path to hdf5.h
HDF5_LIB=
LDFLAGS+=$(HDF5_LIB)

#... Order may matter depending on compiler
LIBS=-llime -fopenmp -lqhg -lDDalphaAMG -lhdf5 -lmxml

# Targets
targets: bench-mg

bench-mg: bench-mg.o mg_interface.o
	$(CC) $(CFLAGS) -o $@ $< mg_interface.o $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o
