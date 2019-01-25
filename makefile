#==============================================================
#
#        MAKEFILE FOR MAGNETISED COLLISIONAL TRANSPORT CODE
#
#        The .o files are compiled separately then linked in
#        the final make MAGNET command.
#
#
#
#===============================================================


include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

CC= mpicc#${PETSC_DIR}/${PETSC_ARCH}/bin/mpicxx
CFLAGS= -fopenmp -O3
LDFLAGS= -fopenmp
OBJS= MAGNET.o Global.o TransportCoeff.o Output.o Initial.o
SRCS= MAGNET.cpp Global.cpp TransportCoeff.cpp Output.cpp Initial.cpp
BINS=MAGNET
INCLUDES= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LIBS=
DIR=~/Desktop/PARAMAGENT
VPATH=src

all: MAGNET

MAGNET: ${OBJS}
	${CC} ${LDFLAGS} ${OBJS} ${PETSC_LIB} -o ${BINS}
	mv MAGNET bin
	rm *.o

MAGNET.o: MAGNET.cpp  Global.h TransportCoeff.h Output.h Initial.h
	${CC} -c ${CFLAGS} $< ${INCLUDES}

Global.o: Global.cpp Global.h
	${CC} -c ${CFLAGS} $<

TransportCoeff.o: TransportCoeff.cpp TransportCoeff.h
	${CC} -c ${CFLAGS} $<

Output.o: Output.cpp Output.h
	${CC} -c ${CFLAGS} $<

Initial.o: Initial.cpp Initial.h
	${CC} -c ${CFLAGS} $<
