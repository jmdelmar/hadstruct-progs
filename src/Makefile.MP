AR=ar
RANLIB=ranlib
CC=mpcc

#... Some module systems will declare this automtacially
HDF5_BASE=/opt/ppc64le/spack/0.9.1/opt/spack/linux-ppc64le/xl-13.1/hdf5-1.8.16-hft56aiahq522onfy5trizdhzkcfjbqq/

#... Adjust according to compiler
CFLAGS=-O2 -I./ -D_GNU_SOURCE 

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
HDF5_INC=-I$(HDF5_BASE)/include
CFLAGS+=$(HDF5_INC)

#... Paths to libraries
LDFLAGS=

#... Path to libqhg.a
QHG_LIB=-L$(HOME)/QuaHoG/lib
LDFLAGS+=$(QHG_LIB)

#... Path to liblime.a
LIME_LIB=-L$(HOME)/loc-install/lib
LDFLAGS+=$(LIME_LIB)

#... Path to libmxml.a
MXML_LIB=-L$(HOME)/loc-install/lib
LDFLAGS+=$(MXML_LIB)

#... Path to libDDalphaAMG.a
MG_LIB=-L$(HOME)/DDalphaAMG/lib
LDFLAGS+=$(MG_LIB)

#... Path to libhdf5.a
HDF5_LIB=-L$(HDF5_BASE)/lib
LDFLAGS+=$(HDF5_LIB)

#... Order may matter depending on compiler
LIBS=-lqhg -lDDalphaAMG -llime -lmxml -lhdf5

# Targets
main: main.o parser.o mg_interface.o
	$(CC) $(CFLAGS) -o main main.o parser.o mg_interface.o $(LDFLAGS) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) *.o
